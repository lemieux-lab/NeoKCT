const RICH_VERSION = 2.0

struct RichKCT{K, Ab, C} <: AbstractKCT{K, Ab}
    kct::NeoKCT{K, Ab, C}
    biotypes::BiotypLayer
    version::Float64
end

RichKCT{K,Ab,C}(kct, biotypes) where {K, Ab<:Alphabet, C} =
    RichKCT{K,Ab,C}(kct, biotypes, RICH_VERSION)

## Delegation ##

Base.length(rich::RichKCT) = length(rich.kct.seqs)
Base.size(rich::RichKCT) = (length(rich),)
idx_prefix_size(rich::RichKCT) = idx_prefix_size(rich.kct)
assemble_count_vector(rich::RichKCT, i::Integer) = assemble_count_vector(rich.kct, i)
compute_index!(rich::RichKCT; kwargs...) = compute_index!(rich.kct; kwargs...)
find_shared_words(rich::RichKCT) = find_shared_words(rich.kct)

Base.findfirst(rich::RichKCT{K,Ab}, key::Kmer{Ab,K}) where {K, Ab} = findfirst(rich.kct, key)
Base.findfirst(rich::RichKCT{K,Ab}, key::UInt64) where {K, Ab} = findfirst(rich.kct, key)

function Base.getindex(rich::RichKCT{K,Ab,C}, i::Integer) where {K, Ab, C}
    kmer, counts = rich.kct[i]
    return kmer => (; counts, biotype = biotype_mask(rich.biotypes, i))
end

Base.getindex(rich::RichKCT{K,Ab,C}, i::UnitRange) where {K,Ab,C} = Tuple(rich[j] for j in i)
Base.getindex(rich::RichKCT{K,Ab,C}, i::V) where {K,Ab,C,V<:AbstractVector} = Tuple(rich[j] for j in i)

# collapse! deduplicates count words without reordering k-mers, so biotype ids are stable.
function collapse!(rich::RichKCT{K,Ab,C}) where {K, Ab, C}
    return RichKCT{K,Ab,C}(collapse!(rich.kct), rich.biotypes)
end

## Biotype merge ##

# Left-join kct k-mers against gidx via a O(n+m) sorted merge walk.
# Matched k-mers get the genomic biotype; unmatched get intergenic (pool[1]).
function add_biotypes(kct::NeoKCT{K,Ab,C}, gidx::GenomicIndex{K,Ab}) where {K, Ab, C}
    n = length(kct.seqs)
    pool = UInt64[INTERGENIC_MASK]
    index_map = Dict{UInt64, UInt16}(INTERGENIC_MASK => UInt16(1))
    ids = fill(UInt16(1), n)

    kct_iter = iterate(kct.seqs)
    gidx_iter = iterate(gidx.seqs)
    i = 1; j = 1

    while !isnothing(kct_iter) && !isnothing(gidx_iter)
        kct_val, kct_state = kct_iter
        gidx_val, gidx_state = gidx_iter
        if kct_val == gidx_val
            ids[i] = _intern_mask!(pool, index_map, get_biotype_mask(gidx, j))
            kct_iter = iterate(kct.seqs, kct_state)
            gidx_iter = iterate(gidx.seqs, gidx_state)
            i += 1; j += 1
        elseif kct_val < gidx_val
            kct_iter = iterate(kct.seqs, kct_state)
            i += 1
        else
            gidx_iter = iterate(gidx.seqs, gidx_state)
            j += 1
        end
    end

    n_intergenic = count(==(UInt16(1)), ids)
    printstyled("Biotype assignment done: $n_intergenic / $n intergenic, $(length(pool)) unique bitmasks\n",
                color=:green)
    return RichKCT{K,Ab,C}(kct, BiotypLayer(ids, pool, gidx.biotype_names))
end
