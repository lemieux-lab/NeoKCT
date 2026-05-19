struct GenomicIndex{K, Ab}
    seqs::DeltaArray
    biotype_ids::Vector{UInt16}  # one per k-mer, index into pool
    pool::Vector{UInt64}
    idx::Pair{Base.RefValue{Int64}, Vector{UnitRange{Int64}}}
    biotype_names::Vector{String}
end

idx_prefix_size(gidx::GenomicIndex) = gidx.idx[1].x

function Base.findfirst(gidx::GenomicIndex{K, Ab}, key::Kmer{Ab, K}) where {K, Ab}
    key_bits = key.data[1]
    idx_key = (key_bits >> (idx_prefix_size(gidx) * bits_per_symbol(Ab()))) + 1
    r = gidx.idx[2][idx_key]
    isempty(r) && return 0
    return searchfirst(gidx.seqs, key_bits, r.start, r.stop)
end

Base.findfirst(gidx::GenomicIndex{K, Ab}, key::UInt64) where {K, Ab} =
    findfirst(gidx, Kmer{Ab, K, 1}(Kmers.unsafe, (key,)))

get_biotype_mask(gidx::GenomicIndex, i::Int) = gidx.pool[gidx.biotype_ids[i]]

function compute_index!(gidx::GenomicIndex{K, Ab}; prefix_size::Int64=4) where {K, Ab}
    start = 1; last_key = 0x0000000000000000
    n = length(gidx.seqs)
    prefix_shift = prefix_size * bits_per_symbol(Ab())
    @showprogress "Computing Genomic Index..." for (i, val) in enumerate(gidx.seqs)
        key = val >> prefix_shift
        if key > last_key
            gidx.idx[2][last_key+1] = start:i
            start = i; last_key = key
        end
    end
    gidx.idx[2][last_key+1] = start:n
    gidx.idx[1].x = prefix_size
end

# Build from a sorted vector of k-mer bits and a parallel vector of biotype bitmasks.
# biotype_names[1] must be "intergenic" (bit 0). Bitmasks are deduplicated internally.
# sorted_kmers must be in ascending order.
function GenomicIndex{K, Ab}(sorted_kmers::Vector{UInt64}, bitmasks::Vector{UInt64},
                              biotype_names::Vector{String};
                              checkpoint_size::Type{<:Unsigned}=UInt64,
                              delta_size::Type{<:Unsigned}=UInt32) where {K, Ab}
    @assert length(sorted_kmers) == length(bitmasks)
    @assert !isempty(biotype_names) && biotype_names[1] == "intergenic"
    pool = UInt64[]
    index_map = Dict{UInt64, UInt16}()
    ids = Vector{UInt16}(undef, length(bitmasks))
    for (i, mask) in enumerate(bitmasks)
        ids[i] = _intern_mask!(pool, index_map, mask)
    end
    seqs = DeltaArray{checkpoint_size, delta_size}(sorted_kmers, DEFAULT_CHECKPOINT_INTERVAL)
    gidx = GenomicIndex{K, Ab}(seqs, ids, pool, Ref(20) => fill(0:-1, 4^15), biotype_names)
    compute_index!(gidx)
    return gidx
end
