const INTERGENIC_MASK = UInt64(1)  # bit 0 reserved for intergenic

struct BiotypLayer
    ids::Vector{UInt16}    # one per k-mer, index into pool (1-based)
    pool::Vector{UInt64}   # deduplicated biotype bitmasks; pool[1] always == INTERGENIC_MASK
    biotype_names::Vector{String}  # biotype_names[b] = name for bit b-1
end

# All-intergenic layer; biotype_names[1] must be "intergenic".
BiotypLayer(n::Int, biotype_names::Vector{String}) =
    BiotypLayer(fill(UInt16(1), n), UInt64[INTERGENIC_MASK], biotype_names)

biotype_mask(l::BiotypLayer, i::Int) = l.pool[l.ids[i]]

function biotype_names_for(l::BiotypLayer, i::Int)
    mask = biotype_mask(l, i)
    return [l.biotype_names[b] for b in eachindex(l.biotype_names) if (mask >> (b - 1)) & 1 == 1]
end

function has_biotype(l::BiotypLayer, i::Int, name::String)
    b = findfirst(==(name), l.biotype_names)
    isnothing(b) && throw(ArgumentError("Unknown biotype: $name"))
    return (biotype_mask(l, i) >> (b - 1)) & 1 == 1
end

# Insert mask into pool if absent; return its 1-based UInt16 index.
function _intern_mask!(pool::Vector{UInt64}, index::Dict{UInt64, UInt16}, mask::UInt64)::UInt16
    return get!(index, mask) do
        push!(pool, mask)
        UInt16(length(pool))
    end
end
