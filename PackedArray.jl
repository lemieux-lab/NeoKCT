"""
    struct PackedArray{T, W <: Unsigned} <: AbstractVector{Vector{T}}
        words::Vector{W}
        bitmap::BitVector
    end

A bitpacked array of type T, using words of size W. Keeps a bitmap to remember element
boundaries within each word. Implements the AbstractVector interface, where each element
is a `Vector{T}` of all values packed into the corresponding word.
"""
struct PackedArray{T, W <: Unsigned} <: AbstractVector{Vector{T}}
    words::Vector{W}
    bitmap::BitVector
end

bitsizeof(e) = sizeof(e) * 8

PackedArray{T, W}() where {T, W<:Unsigned} = PackedArray{T, W}(W[0], BitVector(falses(bitsizeof(W))))
PackedArray{T}() where {T} = PackedArray{T, UInt64}(UInt64[0], BitVector(falses(64)))

## AbstractVector interface ##

Base.size(arr::PackedArray) = (length(arr.words),)
Base.IndexStyle(::Type{<:PackedArray}) = Base.IndexLinear()

# Return all values packed in word i as a Vector{T}
function Base.getindex(arr::PackedArray{T, W}, word_id::Integer) where {T, W<:Unsigned}
    bitmap_word = arr.bitmap[word_bitmap_slice(word_id, W)]
    n = sum(bitmap_word)
    assembly = Vector{T}(undef, n)
    for i in 1:n
        assembly[i] = T(get_val(arr, word_id, i))
    end
    return assembly
end

## Internal helpers ##

# Computes the packed version of elem, returns (packed::W, bitsize::Int)
function prepack(elem::T, ::Type{W}=UInt64) where {T, W<:Unsigned}
    !isbits(elem) && throw("Type $T is not bits type and cannot be packed")
    pack = W(elem & typemax(W)) # /!\ Will round down to typemax, careful
    pack == 0 && return pack, 1
    return pack, bitsizeof(W) - leading_zeros(pack)
end

# Bitmap slice covering the word_id-th word
word_bitmap_slice(word_id::Integer, word_size::DataType=UInt64) =
    (word_id-1)*bitsizeof(word_size)+1 : (word_id-1)*bitsizeof(word_size)+bitsizeof(word_size)

# Last set bit position within the word_id-th word's bitmap slice (0 if none)
function word_last_set(arr::PackedArray{T, W}, word_id::Integer) where {T, W<:Unsigned}
    last_set = findlast(@view arr.bitmap[word_bitmap_slice(word_id, W)])
    return !isnothing(last_set) ? last_set : 0
end

function update_word_bitmap!(arr::PackedArray{T, W}, word_id::Integer, bitsize::Integer) where {T, W<:Unsigned}
    arr.bitmap[word_bitmap_slice(word_id, W)[1]+word_last_set(arr, word_id)+bitsize-1] = true
end

function update_word_bitmap!(arr::PackedArray{T, W}, word_id::Integer, bitsize::Integer, last_set::Integer) where {T, W<:Unsigned}
    arr.bitmap[word_bitmap_slice(word_id, W)[1]+last_set+bitsize-1] = true
end

# Return the elem_id-th value packed in word word_id (as W)
function get_val(arr::PackedArray{T, W}, word_id::Integer, elem_id::Integer) where {T, W<:Unsigned}
    word = arr.words[word_id]
    bitmap_word = arr.bitmap[word_bitmap_slice(word_id, W)]
    sum(bitmap_word) < elem_id && throw("Out of Bound")
    start = 0; finish = 0; counter = 0; ones = 0
    for bit in bitmap_word
        counter += 1
        ones == elem_id && break
        (start, finish, ones) = bit ? (finish, counter, ones+1) : (start, finish, ones)
    end
    return word << start >> (start + (bitsizeof(W)-finish))
end

## Push interface ##

# Pack elem into word at word_id, spilling to a new word if there is no space. Returns arr.
function Base.push!(arr::PackedArray{T1, W}, elem::T2, word_id::Integer) where {T1, T2, W<:Unsigned}
    elem, bitsize = prepack(elem, W)
    last_set = word_last_set(arr, word_id)
    if bitsizeof(W)-last_set >= bitsize
        arr.words[word_id] |= (elem<<(bitsizeof(W)-last_set-bitsize))
        update_word_bitmap!(arr, word_id, bitsize, last_set)
    else
        push!(arr.words, W(0))
        append!(arr.bitmap, falses(bitsizeof(W)))
        push!(arr, elem, length(arr.words))
    end
    return arr
end

# Pack elem into a brand-new word. Returns arr.
function Base.push!(arr::PackedArray{T, W}, elem) where {T, W<:Unsigned}
    elem, _ = prepack(elem, W)
    push!(arr.words, W(0))
    append!(arr.bitmap, falses(bitsizeof(W)))
    push!(arr, elem, length(arr.words))
    return arr
end

# Allocate a new empty word and return its index.
function new_word!(arr::PackedArray{T, W}) where {T, W<:Unsigned}
    push!(arr.words, W(0))
    append!(arr.bitmap, falses(bitsizeof(W)))
    return length(arr.words)
end

## Deduplication ##

# Deduplicates words; returns (deduped_arr, perm, global_perm).
# perm[i] = original word index of the i-th unique word.
# global_perm[i] = index of word i in the deduped array.
function permdedup(arr::PackedArray{T, W}) where {T, W<:Unsigned}
    uniq = Vector{W}()
    indexof = Dict{Pair{W, SubArray}, UInt32}()
    perm = Vector{UInt32}()
    global_perm = Vector{UInt32}(undef, length(arr.words))

    @inbounds @showprogress "Deduping Count Words..." for i in eachindex(arr.words)
        word = arr.words[i]
        bits = @view arr.bitmap[word_bitmap_slice(i, W)]
        j = get(indexof, word=>bits, 0)
        if j == 0
            push!(uniq, word)
            j = length(uniq)
            indexof[word=>bits] = UInt32(j)
            push!(perm, UInt32(i))
        end
        global_perm[i] = j
    end

    bitmap_slices = word_bitmap_slice.(perm, W)
    deduped_arr = PackedArray{T, W}(uniq, arr.bitmap[vcat(bitmap_slices...)])
    return deduped_arr, perm, global_perm
end
