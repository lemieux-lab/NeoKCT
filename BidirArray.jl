
struct BidirArray{T}
    registry::Dict{T, UInt32}
    array::Vector{T}
end

BidirArray(T) = BidirArray(Dict{T, UInt32}(), Vector{T}())
BidirArray{T}() where {T} = BidirArray(T)


function Base.push!(arr::BidirArray{T}, val::T) where {T}
    !isnothing(arr[val]) && return
    push!(arr.array, val)
    arr.registry[val] = UInt32(length(arr.array))
end

Base.getindex(arr::BidirArray{T}, val::T) where {T} = get(arr.registry, val, nothing)
Base.getindex(arr::BidirArray{T}, val::Int) where {T} = get(arr.array, val, nothing)
Base.getindex(arr::BidirArray{T}, vals::V) where {T, V <: AbstractArray{Int}} = getindex.(Ref(arr), vals)
Base.getindex(arr::BidirArray{T}, vals::V) where {T, V <: AbstractRange} = getindex.(Ref(arr), vals)
