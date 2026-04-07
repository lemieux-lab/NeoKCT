if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(".")
end

using Kmers
using BioSequences
using BioSymbols
using ProgressMeter
using Base.Threads
using Dates

include("parallel_sort.jl")
include("PackedArray.jl")
include("JelloFish.jl")

global const VERSION = 1.3

struct NeoKCT{K, Ab<:Alphabet, W<:Unsigned, C}
    seqs::Vector{UInt64}
    n_cids::Vector{UInt16}  # Number of cids per k-mer (essentially delta-encoded offset values of cid positions)
    flat_cids::Vector{UInt32}
    counts::PackedArray{UInt32, W}
    idx::Pair{Base.RefValue{Int64}, Vector{UnitRange{Int64}}}
    samples::Base.RefValue{Int64}
    version::Float64
end

# Convenience constructors
NeoKCT{K,Ab,W,C}(seqs, n_cids, flat_cids, counts, idx, samples) where {K,Ab<:Alphabet,W<:Unsigned,C} =
    NeoKCT{K,Ab,W,C}(seqs, n_cids, flat_cids, counts, idx, samples, VERSION)

NeoKCT{K,Ab,W}() where {K,Ab<:Alphabet,W<:Unsigned} = NeoKCT{K,Ab,W,1}(
    UInt64[], UInt16[], UInt32[],
    PackedArray{UInt32,W}(),
    Ref(20) => fill(0:-1, 4^15),
    Ref(1), VERSION)

# Compute a 1-based start offset into flat_cids for each k-mer from n_cids.
# Returns Vector{UInt64} of length n+1: offsets[i] is the start of k-mer i's cids,
# offsets[n+1] = length(flat_cids) + 1. 
# This is O(n) on n_cids and heavy on memory. Use sparingly.
function _kmer_offsets(n_cids::Vector{UInt16})::Vector{UInt64}
    offsets = Vector{UInt64}(undef, length(n_cids) + 1)
    offsets[1] = 1
    for i in eachindex(n_cids)
        offsets[i+1] = offsets[i] + n_cids[i]
    end
    return offsets
end

# Build KCT of 1 sample from a count hashtable from JelloFish
function NeoKCT{K, Ab, W}(sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet, W<:Unsigned}
    kct = NeoKCT{K, Ab, W}()
    @showprogress desc="Parsing Hash-Table into KCT..." for (k_bits, count) in sample_hashtable
        word_id = push!(kct.counts, UInt64(count))
        push!(kct.seqs, k_bits)
        push!(kct.flat_cids, UInt32(word_id))
        push!(kct.n_cids, UInt16(1))  # each new k-mer starts with exactly 1 cid
    end
    sort!(kct)
    compute_index!(kct)
    return kct
end

idx_prefix_size(kct::NeoKCT) = kct.idx[1].x

# Sorting a KCT is done via psort!, sorts the k-mer sequences, distributing across pivots between threads
function Base.sort!(kct::NeoKCT)
    isempty(kct.seqs) && return kct
    perm = psortperm(kct.seqs)

    kct.seqs[:] = kct.seqs[perm]

    offsets = _kmer_offsets(kct.n_cids)
    new_flat = UInt32[]
    sizehint!(new_flat, length(kct.flat_cids))
    new_n_cids = UInt16[]
    sizehint!(new_n_cids, length(kct.seqs))

    for i in perm
        lo, hi = offsets[i], offsets[i+1] - 1
        append!(new_flat, @view kct.flat_cids[lo:hi])
        push!(new_n_cids, kct.n_cids[i])
    end

    resize!(kct.flat_cids, length(new_flat)); copyto!(kct.flat_cids, new_flat)
    resize!(kct.n_cids, length(new_n_cids)); copyto!(kct.n_cids, new_n_cids)
    return kct
end

# Recomputes the kct.idx based on k-mer prefixes. This speeds up binary search.
function compute_index!(kct::NeoKCT{K, Ab}; prefix_size::Int64=4) where {K, Ab<:Alphabet}
    start = 1
    last_key = 0x0000000000000000
    prefix_shift = prefix_size * bits_per_symbol(Ab())

    @showprogress "Computing Binary Search Index..." for i in eachindex(kct.seqs)
        key = kct.seqs[i] >> prefix_shift
        if key > last_key
            kct.idx[2][last_key+1] = start:i
            start = i
            last_key = key
        end
    end
    kct.idx[2][last_key+1] = start:length(kct.seqs)
    kct.idx[1].x = prefix_size
    return
end

# WARNING: Will render "in place" KCT invalid (needs more work)
function collapse!(kct::NeoKCT{K, Ab, W, C}) where {K, Ab<:Alphabet, W<:Unsigned, C}
    deduped, perms, global_perms = permdedup(kct.counts)
    new_flat = similar(kct.flat_cids)
    @threads for i in eachindex(kct.flat_cids)
        new_flat[i] = UInt32(global_perms[kct.flat_cids[i]])
    end
    printstyled("Collapse done\n", color=:green)
    return NeoKCT{K, Ab, W, C}(kct.seqs, kct.n_cids, new_flat, deduped, kct.idx, kct.samples)
end

# Search for a k-mer in the table. Uses binary search and the k-mer's prefix to search a small region only
function Base.findfirst(kct::NeoKCT{K, Ab}, key::Kmer{Ab, K}) where {K, Ab<:Alphabet}
    idx_key = key.data[1] >> (idx_prefix_size(kct) * bits_per_symbol(Ab())) + 1
    r = kct.idx[2][idx_key]
    i = searchsortedfirst(kct.seqs, key.data[1], r.start, r.stop, Base.Forward)
    return (i in r) && kct.seqs[i] == key.data[1] ? i : 0
end

# From raw UInt64 representing a k-mer stored on a single chunk
# TODO: This should be default, with assembled k-mers calling this given that seqs are directly stored
function Base.findfirst(kct::NeoKCT{K, Ab}, key::UInt64) where {K, Ab<:Alphabet}
    return findfirst(kct, Kmer{Ab, K, 1}(Kmers.unsafe, (key,)))
end

# Push logic for k-mer that was already seen (private, from push!)
function _push_existing_kmer_counts!(kct::NeoKCT{K, Ab, W}, ext_buf::Dict{Int,Vector{UInt32}},
                                      k_pos::Int, shared_words::Set{UInt32}, count::UInt32,
                                      offsets::Vector{UInt64}) where {K, Ab<:Alphabet, W<:Unsigned}
    cur_cids = @view kct.flat_cids[offsets[k_pos] : offsets[k_pos+1]-1]
    last_wid = Int(cur_cids[end])
    total_stored = sum(sum(@view kct.counts.bitmap[word_bitmap_slice(Int(c), W)]) for c in cur_cids)
    missing_counts = kct.samples.x - total_stored
    ext = get!(ext_buf, k_pos, UInt32[])

    if last_wid in shared_words
        if missing_counts > 0
            new_wid = push!(kct.counts, W(0))
            push!(ext, UInt32(new_wid))
            for _ in 2:missing_counts
                next_wid = push!(kct.counts, W(0), new_wid)
                if next_wid != new_wid
                    push!(ext, UInt32(next_wid)); new_wid = next_wid
                end
            end
            word_id = push!(kct.counts, W(count & typemax(W)), new_wid)
            word_id != new_wid && push!(ext, UInt32(word_id))
        else
            word_id = push!(kct.counts, W(count & typemax(W)))
            push!(ext, UInt32(word_id))
        end
    else
        for _ in 1:missing_counts
            new_wid = push!(kct.counts, W(0), last_wid)
            if last_wid != new_wid
                push!(ext, UInt32(new_wid)); last_wid = new_wid
            end
        end
        word_id = push!(kct.counts, W(count & typemax(W)), last_wid)
        last_wid != word_id && push!(ext, UInt32(word_id))
    end
end

# Push logic for k-mer that was never seen (private, from push!)
function _push_new_kmer_counts!(counts::PackedArray{UInt32, W}, prev_samples::Int, count::UInt32) where {W<:Unsigned}
    wid = new_word!(counts)
    chunk_ids = UInt32[wid]
    @inbounds for _ in 1:prev_samples
        new_wid = push!(counts, W(0), wid)
        if new_wid != wid
            push!(chunk_ids, UInt32(new_wid))
            wid = new_wid
        end
    end
    new_wid = push!(counts, W(count & typemax(W)), wid)
    new_wid != wid && push!(chunk_ids, UInt32(new_wid))
    return chunk_ids  # Vector{UInt32} callers append directly to flat_cids
end

# Go through the chunk_ids to find words pointed to by multiple k-mers
# Used to allow word-forking in push! when adding a count to a k-mer but not the other
function find_shared_words(kct::NeoKCT)
    seen = Set{UInt32}()
    shared = Set{UInt32}()
    for cid in kct.flat_cids
        cid in seen ? push!(shared, cid) : push!(seen, cid)
    end
    return shared
end

# Factimily push! function to interface NTuple like a vector for a more "painless" switch to them in chunk_ids
function push(tuple::NTuple{N, T}, val::T) where {N, T}
    return  NTuple{N+1, T}((tuple..., val))
end

# Main push! from a sample to an existing KCT
function Base.push!(kct::NeoKCT{K, Ab, W}, sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet, W<:Unsigned}
    shared_words = find_shared_words(kct)
    offsets = _kmer_offsets(kct.n_cids)  # compute once for all per-kmer lookups
    ext_buf = Dict{Int, Vector{UInt32}}()
    new_seqs = UInt64[]
    new_cids = Vector{UInt32}[]

    @showprogress desc="Adding Sample $(kct.samples.x+1) to Table..." for (k_bits, count) in sample_hashtable
        tmp_seq = Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,))
        k_pos = findfirst(kct, tmp_seq)

        if k_pos == 0
            push!(new_seqs, k_bits)
            push!(new_cids, _push_new_kmer_counts!(kct.counts, kct.samples.x, count))
        else
            _push_existing_kmer_counts!(kct, ext_buf, k_pos, shared_words, count, offsets)
        end
    end

    _merge_and_sort!(kct, ext_buf, new_seqs, new_cids, offsets)
    compute_index!(kct)
    kct.samples.x += 1
end

function _merge_and_sort!(kct::NeoKCT, ext_buf::Dict{Int,Vector{UInt32}},
                           new_seqs::Vector{UInt64}, new_cids::Vector{Vector{UInt32}},
                           offsets::Vector{UInt64})
    n_existing = length(kct.seqs)
    n_new = length(new_seqs)
    n_total = n_existing + n_new

    all_seqs = vcat(kct.seqs, new_seqs)
    perm = psortperm(all_seqs)

    new_flat = UInt32[]
    sizehint!(new_flat, length(kct.flat_cids) + sum(length, values(ext_buf); init=0) +
                        sum(length, new_cids; init=0))
    new_n_cids = UInt16[]
    sizehint!(new_n_cids, n_total)

    for i in perm
        cids_start = length(new_flat)
        if i <= n_existing
            append!(new_flat, @view kct.flat_cids[offsets[i] : offsets[i+1]-1])
            haskey(ext_buf, i) && append!(new_flat, ext_buf[i])
        else
            append!(new_flat, new_cids[i - n_existing])
        end
        push!(new_n_cids, UInt16(length(new_flat) - cids_start))  # always small, no overflow
    end

    resize!(kct.seqs, n_total); copyto!(kct.seqs, all_seqs[perm])
    resize!(kct.flat_cids, length(new_flat)); copyto!(kct.flat_cids, new_flat)
    resize!(kct.n_cids, length(new_n_cids)); copyto!(kct.n_cids, new_n_cids)
end

function assemble_count_vector(kct::NeoKCT{K, Ab, W}, i::Integer) where {K, Ab<:Alphabet, W<:Unsigned}
    lo = 1 + Int64(sum(@view kct.n_cids[1:i-1]))
    cids = @view kct.flat_cids[lo : lo + kct.n_cids[i] - 1]
    counts = reduce(vcat, [assemble_word(kct.counts, Int(c)) for c in cids])
    n_missing = kct.samples.x - length(counts)
    return n_missing > 0 ? vcat(counts, zeros(W, n_missing)) : counts
end

function Base.getindex(kct::NeoKCT{K, Ab, W}, i::Integer) where {K, Ab<:Alphabet, W<:Unsigned}
    seq = Kmer{Ab, K, 1}(Kmers.unsafe, (kct.seqs[i],))
    return seq => assemble_count_vector(kct, i)
end

function build_kct(sample::String, K::Int=30, chunks::Int = 500_000; word_size::DataType=UInt128, collapse::Bool=true)
    kct = NeoKCT{K÷3, AAAlphabet, word_size}(jello_superthreaded_hash(sample, K, chunks))
    kct = collapse ? collapse!(kct) : kct
    return kct
end

function build_kct(samples::AbstractVector{String}, K::Int=30, chunks::Int = 500_000; word_size::DataType=UInt128, save_at_samples::AbstractVector{Int}=Int[], save_path::String = "", collapse_every::Int=1, benchmark_every::Int=5, full_pointer_walkthrough::Bool=false)
    kct = NeoKCT{K÷3, AAAlphabet, word_size}(jello_superthreaded_hash(popfirst!(samples), K, chunks))
    mkpath(save_path*"benchmarks/")
    for (i, sample) in enumerate(samples)
        sample_hashtable = jello_superthreaded_hash(sample, K, chunks)
        push!(kct, sample_hashtable)
        kct = ((i+1) % collapse_every == 0) ? collapse!(kct) : kct 
        i+1 in save_at_samples && write_kct(kct, save_path*"$(day(now()))_$(month(now()))_$(year(now()))_[$(i+1)_samples]_neokct.kct")
        ((i+1) % benchmark_every == 0) && benchmark_kct(kct, save_path*"benchmarks/", full_pointer_walkthrough=full_pointer_walkthrough)
    end
    return kct
end

include("KCTBenchmarker.jl")
include("KCTLoader.jl")

### EXPERIMENTAL ### 

# temporary path to quickly test on all samples
samples = readlines(open("/u/jacquinn/phd_stuff/data/tcga_fastqs.paths", "r"))
