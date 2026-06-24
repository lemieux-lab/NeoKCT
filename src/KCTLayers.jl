using Kmers
using BioSequences
using BioSymbols
using ProgressMeter
using Base.Threads
using NArrays
using Dates

include("JelloFish.jl")

abstract type AbstractLayer end

## K-mer Layer ## 

struct KmerLayer{K, Ab <: Alphabet, C <: Unsigned, D <: Unsigned} <: AbstractLayer
    seqs::DeltaArray{C, D}
    idx::Pair{Base.RefValue{Int64}, Vector{UnitRange{Int64}}}
end

idx_prefix_size(klayer::KmerLayer) = klayer.idx[1].x

function Base.getindex(klayer::KmerLayer, i::UnitRange)
    return Tuple(klayer[j] for j in i)
end

function Base.getindex(klayer::KmerLayer, i::AbstractVector)
    return Tuple(klayer[j] for j in i)
end

const DEFAULT_CHECKPOINT_INTERVAL = 256

Base.length(kl::KmerLayer) = length(kl.seqs)

function Base.getindex(kl::KmerLayer{K, Ab}, i::Integer) where {K, Ab<:Alphabet}
    return Kmer{Ab, K, 1}(Kmers.unsafe, (kl.seqs[i],))
end

KmerLayer{K, Ab}(; checkpoint_size::Type{C}=UInt64,
                   delta_size::Type{D}=UInt32) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KmerLayer{K, Ab, C, D}(DeltaArray{checkpoint_size, delta_size}(DEFAULT_CHECKPOINT_INTERVAL),
                            Ref(20) => fill(0:-1, 4^15))

function Base.findfirst(kl::KmerLayer{K, Ab}, key::Kmer{Ab, K}) where {K, Ab<:Alphabet}
    key_bits = key.data[1]
    idx_key = (key_bits >> (idx_prefix_size(kl) * bits_per_symbol(Ab()))) + 1
    r = kl.idx[2][idx_key]
    isempty(r) && return 0
    return searchfirst(kl.seqs, key_bits, r.start, r.stop)
end

Base.findfirst(kl::KmerLayer{K, Ab}, key::UInt64) where {K, Ab<:Alphabet} =
    findfirst(kl, Kmer{Ab, K, 1}(Kmers.unsafe, (key,)))

function compute_index!(kl::KmerLayer{K, Ab}; prefix_size::Int64=4) where {K, Ab<:Alphabet}
    start = 1
    last_key = 0x0000000000000000
    n = length(kl.seqs)
    prefix_shift = prefix_size * bits_per_symbol(Ab())
    @showprogress "Computing Binary Search Index..." for (i, val) in enumerate(kl.seqs)
        key = val >> prefix_shift
        if key > last_key
            kl.idx[2][last_key+1] = start:i
            start = i
            last_key = key
        end
    end
    kl.idx[2][last_key+1] = start:n
    kl.idx[1].x = prefix_size
end

## Counts Layer ##

struct CountsLayer <: AbstractLayer
    n_cids::Vector{UInt16}  # number of cids per k-mer (cumsum to reconstruct offset vector)
    flat_cids::Vector{UInt32}
    counts::PackedArray{UInt32}
    samples::Base.RefValue{Int64}
end

Base.length(cl::CountsLayer) = length(cl.n_cids)

CountsLayer(; word_size::Type{<:Unsigned}=UInt128) =
    CountsLayer(UInt16[], UInt32[], PackedArray{UInt32, word_size}(), Ref(0))

function assemble_count_vector(cl::CountsLayer, i::Integer)
    lo = 1 + Int64(sum(@view cl.n_cids[1:i-1]))
    cids = @view cl.flat_cids[lo : lo + cl.n_cids[i] - 1]
    counts = reduce(vcat, [cl.counts[Int(c)] for c in cids])
    n_missing = cl.samples.x - length(counts)
    return n_missing > 0 ? vcat(counts, zeros(eltype(cl.counts.words), n_missing)) : counts
end

Base.getindex(cl::CountsLayer, i::Integer) = assemble_count_vector(cl, i)

function find_shared_words(cl::CountsLayer)
    seen = Set{UInt32}()
    shared = Set{UInt32}()
    for cid in cl.flat_cids
        cid in seen ? push!(shared, cid) : push!(seen, cid)
    end
    return shared
end

function _kmer_offsets(n_cids::Vector{UInt16})::Vector{UInt64}
    offsets = Vector{UInt64}(undef, length(n_cids) + 1)
    offsets[1] = 1
    for i in eachindex(n_cids)
        offsets[i+1] = offsets[i] + n_cids[i]
    end
    return offsets
end

function collapse!(cl::CountsLayer)
    deduped, _, global_perms = permdedup(cl.counts)
    new_flat = similar(cl.flat_cids)
    @threads for i in eachindex(cl.flat_cids)
        new_flat[i] = UInt32(global_perms[cl.flat_cids[i]])
    end
    printstyled("Collapse done\n", color=:green)
    return CountsLayer(cl.n_cids, new_flat, deduped, cl.samples)
end

## Biotype Layer ##

const INTERGENIC_MASK = UInt64(1)  # bit 0 reserved for intergenic

struct BiotypLayer <: AbstractLayer
    ids::Vector{UInt16}    # one per k-mer, index into pool (1-based)
    pool::Vector{UInt64}   # deduplicated biotype bitmasks; pool[1] always == INTERGENIC_MASK
    biotype_names::Vector{String}  # biotype_names[b] = name for bit b-1
end

# All-intergenic layer; biotype_names[1] must be "intergenic".
BiotypLayer(n::Int, biotype_names::Vector{String}) =
    BiotypLayer(fill(UInt16(1), n), UInt64[INTERGENIC_MASK], biotype_names)

Base.length(bl::BiotypLayer) = length(bl.ids)

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

## KCT Wrapper ##

struct KCT{K, Ab <: Alphabet, Counts <: Union{CountsLayer, Nothing}, Biotype <: Union{BiotypLayer, Nothing}, C <: Unsigned, D <: Unsigned}
    kmer::KmerLayer{K, Ab, C, D}
    counts::Counts
    biotype::Biotype
end

KCT(kl::KmerLayer{K, Ab, C, D}) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, Nothing, Nothing, C, D}(kl, nothing, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, Nothing, C, D}(kl, cl, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, BiotypLayer, C, D}(kl, cl, bl)
KCT(kl::KmerLayer{K, Ab, C, D}, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, Nothing, BiotypLayer, C, D}(kl, nothing, bl)

## KCT Methods ##

Base.length(kct::KCT) = length(kct.kmer)
Base.size(kct::KCT) = (length(kct),)
idx_prefix_size(kct::KCT) = idx_prefix_size(kct.kmer)
compute_index!(kct::KCT; kwargs...) = compute_index!(kct.kmer; kwargs...)

Base.findfirst(kct::KCT{K, Ab}, key::Kmer{Ab, K}) where {K, Ab} = findfirst(kct.kmer, key)
Base.findfirst(kct::KCT{K, Ab}, key::UInt64) where {K, Ab} = findfirst(kct.kmer, key)

function Base.getindex(kct::KCT{K, Ab, Nothing, Nothing}, i::Integer) where {K, Ab}
    return kct.kmer[i]
end

function Base.getindex(kct::KCT{K, Ab, CountsLayer, Nothing}, i::Integer) where {K, Ab}
    return kct.kmer[i] => kct.counts[i]
end

function Base.getindex(kct::KCT{K, Ab, CountsLayer, BiotypLayer}, i::Integer) where {K, Ab}
    return kct.kmer[i] => (; counts=kct.counts[i], biotype=biotype_mask(kct.biotype, i))
end

function Base.getindex(kct::KCT{K, Ab, Nothing, BiotypLayer}, i::Integer) where {K, Ab}
    return kct.kmer[i] => biotype_mask(kct.biotype, i)
end

Base.getindex(kct::KCT, i::UnitRange) = Tuple(kct[j] for j in i)
Base.getindex(kct::KCT, i::AbstractVector) = Tuple(kct[j] for j in i)

function collapse!(kct::KCT{K, Ab, CountsLayer, Nothing}) where {K, Ab}
    return KCT(kct.kmer, collapse!(kct.counts))
end

function collapse!(kct::KCT{K, Ab, CountsLayer, BiotypLayer}) where {K, Ab}
    return KCT(kct.kmer, collapse!(kct.counts), kct.biotype)
end

## KCT Sample Building ##

function KCT{K, Ab}(sample_hashtable::Dict{UInt64, UInt32};
                    checkpoint_size::Type{C}=UInt64,
                    delta_size::Type{D}=UInt32,
                    word_size::Type{<:Unsigned}=UInt128) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned}
    kl = KmerLayer{K, Ab}(checkpoint_size=checkpoint_size, delta_size=delta_size)
    cl = CountsLayer(word_size=word_size)
    tmp_seqs = Vector{UInt64}(undef, length(sample_hashtable))
    sizehint!(cl.flat_cids, length(sample_hashtable))
    sizehint!(cl.n_cids, length(sample_hashtable))
    @showprogress desc="Parsing Hash-Table into KCT..." for (j, (k_bits, count)) in enumerate(sample_hashtable)
        push!(cl.counts, UInt64(count))
        tmp_seqs[j] = k_bits
        push!(cl.flat_cids, UInt32(lastindex(cl.counts)))
        push!(cl.n_cids, UInt16(1))
    end
    perm = psortperm(tmp_seqs)
    new_flat = similar(cl.flat_cids); new_n_cids = similar(cl.n_cids)
    for (j, i) in enumerate(perm)
        new_flat[j]   = cl.flat_cids[i]
        new_n_cids[j] = cl.n_cids[i]
    end
    copyto!(cl.flat_cids, new_flat); copyto!(cl.n_cids, new_n_cids)
    encode!(kl.seqs, tmp_seqs[perm])
    cl.samples.x = 1
    kct = KCT(kl, cl)
    compute_index!(kct.kmer)
    return kct
end

## KCT push! helpers ##

function _push_new_kmer_counts!(cl::CountsLayer, count::UInt32)
    W = eltype(cl.counts.words)
    wid = new_word!(cl.counts)
    chunk_ids = UInt32[wid]
    @inbounds for _ in 1:cl.samples.x
        push!(cl.counts, W(0), wid)
        if lastindex(cl.counts) != wid
            wid = lastindex(cl.counts)
            push!(chunk_ids, UInt32(wid))
        end
    end
    push!(cl.counts, W(count & typemax(W)), wid)
    lastindex(cl.counts) != wid && push!(chunk_ids, UInt32(lastindex(cl.counts)))
    return chunk_ids
end

function _push_existing_kmer_counts!(cl::CountsLayer, ext_buf::Dict{Int, Vector{UInt32}},
                                      k_pos::Int, shared_words::Set{UInt32}, count::UInt32,
                                      offsets::Vector{UInt64})
    W = eltype(cl.counts.words)
    cur_cids = @view cl.flat_cids[offsets[k_pos] : offsets[k_pos+1]-1]
    last_wid = Int(cur_cids[end])
    total_stored = sum(length(cl.counts[Int(c)]) for c in cur_cids)
    missing_counts = cl.samples.x - total_stored
    ext = get!(ext_buf, k_pos, UInt32[])

    if last_wid in shared_words
        if missing_counts > 0
            push!(cl.counts, W(0))
            new_wid = lastindex(cl.counts)
            push!(ext, UInt32(new_wid))
            for _ in 2:missing_counts
                push!(cl.counts, W(0), new_wid)
                next_wid = lastindex(cl.counts)
                if next_wid != new_wid
                    push!(ext, UInt32(next_wid)); new_wid = next_wid
                end
            end
            push!(cl.counts, W(count & typemax(W)), new_wid)
            word_id = lastindex(cl.counts)
            word_id != new_wid && push!(ext, UInt32(word_id))
        else
            push!(cl.counts, W(count & typemax(W)))
            push!(ext, UInt32(lastindex(cl.counts)))
        end
    else
        for _ in 1:missing_counts
            push!(cl.counts, W(0), last_wid)
            new_wid = lastindex(cl.counts)
            if last_wid != new_wid
                push!(ext, UInt32(new_wid)); last_wid = new_wid
            end
        end
        push!(cl.counts, W(count & typemax(W)), last_wid)
        word_id = lastindex(cl.counts)
        last_wid != word_id && push!(ext, UInt32(word_id))
    end
end

function _merge_and_sort!(kl::KmerLayer, cl::CountsLayer, bl::Union{BiotypLayer, Nothing},
                           ext_buf::Dict{Int, Vector{UInt32}}, new_seqs::Vector{UInt64},
                           new_cids::Vector{Vector{UInt32}}, offsets::Vector{UInt64})
    old_seqs = collect(kl.seqs)
    n_existing = length(old_seqs)
    all_seqs = vcat(old_seqs, new_seqs)
    perm = psortperm(all_seqs)

    new_flat = UInt32[]
    sizehint!(new_flat, length(cl.flat_cids) + sum(length, values(ext_buf); init=0) +
                        sum(length, new_cids; init=0))
    new_n_cids = UInt16[]
    sizehint!(new_n_cids, length(all_seqs))
    new_ids = isnothing(bl) ? nothing : Vector{UInt16}(undef, length(all_seqs))

    for (j, i) in enumerate(perm)
        cids_start = length(new_flat)
        if i <= n_existing
            append!(new_flat, @view cl.flat_cids[offsets[i] : offsets[i+1]-1])
            haskey(ext_buf, i) && append!(new_flat, ext_buf[i])
            !isnothing(new_ids) && (new_ids[j] = bl.ids[i])
        else
            append!(new_flat, new_cids[i - n_existing])
            !isnothing(new_ids) && (new_ids[j] = UInt16(1))  # intergenic for new k-mers
        end
        push!(new_n_cids, UInt16(length(new_flat) - cids_start))
    end

    if !isnothing(bl)
        resize!(bl.ids, length(all_seqs))
        copyto!(bl.ids, new_ids)
    end
    encode!(kl.seqs, all_seqs[perm])
    resize!(cl.flat_cids, length(new_flat)); copyto!(cl.flat_cids, new_flat)
    resize!(cl.n_cids, length(new_n_cids)); copyto!(cl.n_cids, new_n_cids)
end

## KCT push! and sort! ##

function Base.push!(kct::KCT{K, Ab, CountsLayer, Biotype},
                    sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet, Biotype}
    cl = kct.counts
    shared_words = find_shared_words(cl)
    offsets = _kmer_offsets(cl.n_cids)
    ext_buf = Dict{Int, Vector{UInt32}}()
    new_seqs = UInt64[]
    new_cids = Vector{UInt32}[]

    @showprogress desc="Adding Sample $(cl.samples.x + 1) to Table..." for (k_bits, count) in sample_hashtable
        k_pos = findfirst(kct.kmer, Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,)))
        if k_pos == 0
            push!(new_seqs, k_bits)
            push!(new_cids, _push_new_kmer_counts!(cl, count))
        else
            _push_existing_kmer_counts!(cl, ext_buf, k_pos, shared_words, count, offsets)
        end
    end

    _merge_and_sort!(kct.kmer, cl, kct.biotype, ext_buf, new_seqs, new_cids, offsets)
    compute_index!(kct.kmer)
    cl.samples.x += 1
end

function Base.sort!(kct::KCT{K, Ab, CountsLayer, Biotype}) where {K, Ab, Biotype}
    isempty(kct.kmer.seqs) && return kct
    seqs = collect(kct.kmer.seqs)
    perm = psortperm(seqs)
    encode!(kct.kmer.seqs, seqs[perm])

    cl = kct.counts
    offsets = _kmer_offsets(cl.n_cids)
    new_flat = UInt32[]; sizehint!(new_flat, length(cl.flat_cids))
    new_n_cids = UInt16[]; sizehint!(new_n_cids, length(cl.n_cids))
    for i in perm
        lo, hi = offsets[i], offsets[i+1] - 1
        append!(new_flat, @view cl.flat_cids[lo:hi])
        push!(new_n_cids, cl.n_cids[i])
    end
    resize!(cl.flat_cids, length(new_flat)); copyto!(cl.flat_cids, new_flat)
    resize!(cl.n_cids, length(new_n_cids)); copyto!(cl.n_cids, new_n_cids)

    if !isnothing(kct.biotype)
        new_ids = kct.biotype.ids[perm]
        copyto!(kct.biotype.ids, new_ids)
    end
    return kct
end

## KCT Genomic Index ##

# Build a KCT{K, Ab, Nothing, BiotypLayer} from a sorted k-mer vector and a parallel bitmask vector.
# Replaces the old GenomicIndex{K,Ab} constructor.
function KCT{K, Ab}(sorted_kmers::Vector{UInt64}, bitmasks::Vector{UInt64},
                    biotype_names::Vector{String};
                    checkpoint_size::Type{C}=UInt64,
                    delta_size::Type{D}=UInt32) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned}
    @assert length(sorted_kmers) == length(bitmasks)
    @assert !isempty(biotype_names) && biotype_names[1] == "intergenic"
    pool = UInt64[]
    index_map = Dict{UInt64, UInt16}()
    ids = Vector{UInt16}(undef, length(bitmasks))
    for (i, mask) in enumerate(bitmasks)
        ids[i] = _intern_mask!(pool, index_map, mask)
    end
    kl = KmerLayer{K, Ab, C, D}(DeltaArray{checkpoint_size, delta_size}(sorted_kmers, DEFAULT_CHECKPOINT_INTERVAL),
                                  Ref(20) => fill(0:-1, 4^15))
    kct = KCT(kl, BiotypLayer(ids, pool, biotype_names))
    compute_index!(kct.kmer)
    return kct
end

# Left-join kct k-mers against gidx via an O(n+m) sorted merge walk.
# Matched k-mers get the genomic biotype; unmatched k-mers stay intergenic.
function add_biotypes(kct::KCT{K, Ab, CountsLayer, Nothing, C, D},
                      gidx::KCT{K, Ab, Nothing, BiotypLayer}) where {K, Ab, C<:Unsigned, D<:Unsigned}
    n = length(kct.kmer)
    pool = UInt64[INTERGENIC_MASK]
    index_map = Dict{UInt64, UInt16}(INTERGENIC_MASK => UInt16(1))
    ids = fill(UInt16(1), n)

    kct_iter = iterate(kct.kmer.seqs)
    gidx_iter = iterate(gidx.kmer.seqs)
    i = 1; j = 1

    while !isnothing(kct_iter) && !isnothing(gidx_iter)
        kct_val, kct_state = kct_iter
        gidx_val, gidx_state = gidx_iter
        if kct_val == gidx_val
            ids[i] = _intern_mask!(pool, index_map, biotype_mask(gidx.biotype, j))
            kct_iter = iterate(kct.kmer.seqs, kct_state)
            gidx_iter = iterate(gidx.kmer.seqs, gidx_state)
            i += 1; j += 1
        elseif kct_val < gidx_val
            kct_iter = iterate(kct.kmer.seqs, kct_state)
            i += 1
        else
            gidx_iter = iterate(gidx.kmer.seqs, gidx_state)
            j += 1
        end
    end

    n_intergenic = count(==(UInt16(1)), ids)
    printstyled("Biotype assignment done: $n_intergenic / $n intergenic, $(length(pool)) unique bitmasks\n",
                color=:green)
    return KCT(kct.kmer, kct.counts, BiotypLayer(ids, pool, gidx.biotype.biotype_names))
end

include("GenomicIndexBuilder.jl")
include("KCTBenchmarker.jl")

function build_kct(samples::AbstractVector{String}, K::Int=30, chunks::Int=500_000;
                   word_size::Type{<:Unsigned}=UInt128, checkpoint_size::Type{<:Unsigned}=UInt64,
                   delta_size::Type{<:Unsigned}=UInt32, save_at_samples::AbstractVector{Int}=Int[],
                   save_path::String="", collapse_every::Int=1, benchmark_every::Int=5,
                   full_pointer_walkthrough::Bool=false)
    kct = KCT{K÷3, AAAlphabet}(jello_superthreaded_hash(popfirst!(samples), K, chunks),
                                checkpoint_size=checkpoint_size,
                                delta_size=delta_size,
                                word_size=word_size)
    length(samples) != 0 && mkpath(save_path * "benchmarks/")
    for (i, sample) in enumerate(samples)
        push!(kct, jello_superthreaded_hash(sample, K, chunks))
        kct = ((i + 1) % collapse_every == 0) ? collapse!(kct) : kct
        (i + 1) in save_at_samples && write_kct(kct, save_path *
            "$(day(now()))_$(month(now()))_$(year(now()))_[$(i+1)_samples]_neokct.kct")
        ((i + 1) % benchmark_every == 0) && benchmark_kct(kct, save_path * "benchmarks/",
            full_pointer_walkthrough=full_pointer_walkthrough)
    end
    return kct
end

include("KCTLoader.jl")

# temporary path to quickly test on all samples
samples = readlines(open("/u/jacquinn/phd_stuff/data/tcga_fastqs.paths", "r"));