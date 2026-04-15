include("parallel_sort.jl")

"""
    struct DeltaArrayPackedArray{C, D} checkpoints::Vector{C} deltas::Vector{D} regular_cp_idx::Vector{UInt64} checkpoint_interval::Int end -> Return type

A delta-encoded array of C, with deltas kept in type D. Only useful if sizeof(D) < sizeof(C).
Checkpoints reduce search time and are found at "checkpoint_interval". Stored as 0s on deltas.
If (element n+1 - element n) > typemax(D), element n+1 becomes a checkpoint instead (indexing kept with regular_cp_index).
"""
struct DeltaArray{C, D} <: AbstractVector{C}
    checkpoints::Vector{C}
    deltas::Vector{D}
    regular_cp_idx::Vector{UInt64}
    checkpoint_interval::Int
end

# Empty constructor
DeltaArray{C, D}(checkpoint_interval::Int=256) where {C, D} =
    DeltaArray(C[], D[], UInt64[], checkpoint_interval)

DeltaArray(checkpoint_interval::Int=256) =
    DeltaArray{UInt64, UInt32}(checkpoint_interval)

# Build from a sorted Vector
DeltaArray{C, D}(arr::Vector{C}, checkpoint_interval::Int=256) where {C, D} =
    DeltaArray(_delta_encode(arr, checkpoint_interval, D)..., checkpoint_interval)

DeltaArray(arr::Vector{C}, checkpoint_interval::Int=256) where {C} =
    DeltaArray(_delta_encode(arr, checkpoint_interval)..., checkpoint_interval)

# Unsorted helper
DeltaArray(arr::Vector{C}, checkpoint_interval::Int=256; presorted::Bool=false) where {C} = 
    DeltaArray(_delta_encode(presorted ? arr : psort!(arr), checkpoint_interval)..., checkpoint_interval)

## Encoding / decoding ##

function _delta_encode(arr::Vector{C}, inter::Int, ::Type{D}=UInt32) where {C, D}
    n = length(arr)
    n == 0 && return C[], D[], UInt64[]
    checkpoints = C[]
    deltas = zeros(D, n)
    regular_cp_idx = Vector{UInt64}(undef, cld(n, inter))
    for i in eachindex(arr)
        is_regular = (i - 1) % inter == 0
        overflow = i > 1 && (arr[i] - arr[i-1]) > typemax(D)
        if is_regular || overflow
            is_regular && (regular_cp_idx[(i-1) ÷ inter + 1] = length(checkpoints) + 1)
            push!(checkpoints, arr[i])
            deltas[i] = D(0)  # sentinel
        else
            deltas[i] = D(arr[i] - arr[i-1])
        end
    end
    return checkpoints, deltas, regular_cp_idx
end

# Re-encode in-place from a new sorted flat vector (used after sort/merge).
function encode!(a::DeltaArray{C, D}, arr::Vector{C}) where {C, D}
    cp, dl, rci = _delta_encode(arr, a.checkpoint_interval, D)
    resize!(a.checkpoints, length(cp)); copyto!(a.checkpoints, cp)
    resize!(a.deltas, length(dl)); copyto!(a.deltas, dl)
    resize!(a.regular_cp_idx, length(rci)); copyto!(a.regular_cp_idx, rci)
end

## AbstractVector Interface ##

Base.length(a::DeltaArray) = length(a.deltas)
Base.size(a::DeltaArray) = (length(a.deltas),)
Base.IndexStyle(::Type{<:DeltaArray}) = Base.IndexLinear()

# Random access: O(checkpoint_interval) worst case.
function Base.getindex(a::DeltaArray{C}, i::Integer)::C where {C}
    @boundscheck checkbounds(a, i)
    inter = a.checkpoint_interval
    k = (i - 1) ÷ inter
    cp = a.regular_cp_idx[k + 1]
    v = a.checkpoints[cp]
    for j in (k * inter + 2):i
        if a.deltas[j] == 0
            cp += 1
            v = a.checkpoints[cp]
        else
            v += a.deltas[j]
        end
    end
    return v
end

# Sequential iteration: O(1) amortised. Avoids re-decoding from checkpoint each step.
# State = (next_index, checkpoint_index, current_value).
function Base.iterate(a::DeltaArray{C}, state=(1, 0, C(0))) where {C}
    i, cp, v = state
    i > length(a) && return nothing
    if a.deltas[i] == 0
        cp += 1
        v = a.checkpoints[cp]
    else
        v += a.deltas[i]
    end
    return v, (i + 1, cp, v)
end

## Sorted search ##

# Return the first index in [lo, hi] whose decoded value equals target,
# or 0 if not found (exits early if decoded value exceeds target).
# Starts from the regular checkpoint covering lo. O(C + hi-lo) total.
function searchfirst(a::DeltaArray{C}, target::C, lo::Int, hi::Int)::Int where {C}
    isempty(a) && return 0
    inter = a.checkpoint_interval
    k = (lo - 1) ÷ inter
    cp = a.regular_cp_idx[k + 1]
    v = a.checkpoints[cp]
    chunk_start = k * inter + 1
    for j in chunk_start:hi
        if j > chunk_start
            if a.deltas[j] == 0
                cp += 1
                v = a.checkpoints[cp]
            else
                v  += a.deltas[j]
            end
        end
        j < lo && continue
        v == target && return j
        v >  target && return 0
    end
    return 0
end