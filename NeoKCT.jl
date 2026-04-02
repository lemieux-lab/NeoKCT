using Pkg
Pkg.activate(".")

using Kmers
using BioSequences
using JSON
using ProgressMeter
using Base.Threads
using Dates
using GZip
using CairoMakie

include("parallel_sort.jl")
include("PackedArray.jl")
include("AAAlphabet.jl")
include("BidirArray.jl")
include("JelloFish.jl")

global const VERSION = 1.1

"""
An element in a K-mer Count Table that keeps a K-mer Sequence of a given Alphabet, and indexes towards chunks
of counts or symbols.

   K_Element{K, Ab<:Alphabet}(seq::Kmer, chunk_ids::Vector{UInt32})
"""
struct K_Element{K, Ab<:Alphabet, C}
    seq::Kmer{Ab, K, C}
    chunk_ids::NTuple
end

K_Element{K, Ab}(seq::S) where {K, Ab<:Alphabet, S<:Kmer{Ab, K}} = K_Element{K, Ab, 1}(seq, NTuple{0, UInt32}())
K_Element{K, Ab}(seq::S) where {K, Ab<:Alphabet, S<:AbstractString} = K_Element{K, Ab, 1}(Kmer{Ab, K, 1}(seq), NTuple{0, UInt32}())
K_Element{K, Ab}(seq::S, chunk_ids::NTuple) where {K, Ab<:Alphabet, S<:Kmer{Ab, K}} = K_Element{K, Ab, 1}(seq, chunk_ids)

"""
A sorted k-mer count table implementing bitpacked count arrays.

    NeoKCT{K, Ab<:Alphabet, W<:Unsigned}()
"""
struct NeoKCT{K, Ab<:Alphabet, W<:Unsigned, C}
    table::Vector{K_Element{K, Ab, C}}
    counts::PackedArray{UInt32, W}
    idx::Pair{Base.RefValue{Int64} , Vector{UnitRange{Int64}}}
    samples::Base.RefValue{Int64}  # Send me to jail for this one
    version::Float64
end

NeoKCT{K, Ab, W, C}(table, counts, idx, samples) where {K, Ab<:Alphabet, W<:Unsigned, C} = NeoKCT(table, counts, idx, samples, VERSION)
NeoKCT{K, Ab, W}() where {K, Ab<:Alphabet, W<:Unsigned} = NeoKCT(K_Element{K, Ab, 1}[], PackedArray{UInt32, W}(), Ref(20)=>fill(0:-1, 4^15), Ref(1), VERSION)


# Comparison rules between types that can be sorted here
Base.isless(a::K_Element, b::K_Element) = a.seq < b.seq 
Base.isless(a::Kmer, b::K_Element) = a < b.seq
Base.isless(a::K_Element, b::Kmer) = a.seq < b
Base.:(==)(a::K_Element, b::K_Element) = a.seq == b.seq

Base.isless(a::K_Element, b::UInt) = a.seq.data[1] < b 
Base.isless(a::UInt, b::K_Element) = a < b.seq.data[1]
Base.:(==)(a::K_Element, b::UInt) = a.seq.data[1] == b
Base.:(==)(a::UInt, b::K_Element) = a == b.seq.data[1]

# Sorting a kct always calls the parallelised sorter on the table
Base.sort!(kct::NeoKCT) = psort!(kct.table)

# Recomputes the kct.idx based on k-mer prefixes. This speeds up binary search.
function compute_index!(kct::NeoKCT{K, Ab}; prefix_size::Int64=4) where {K, Ab<:Alphabet}
    start = 1
    last_key = 0x0000000000000000
    symbol_size = bits_per_symbol(Ab())
    prefix_shift = prefix_size*symbol_size
    
    @showprogress "Computing Binary Search Index..." for i in eachindex(kct.table)
        key = kct.table[i].seq.data[1] >> prefix_shift
        if key > last_key
            # println(last_key)
            # last_key+1 > length(kct.idx) && append!(kct.idx, fill(0:-1, length(kct.idx)))  # Currently memory leaks
            kct.idx[2][last_key+1] = start:i
            start = i
            last_key = key
        end
    end
    kct.idx[2][last_key+1] = start:length(kct.table)
    kct.idx[1].x = prefix_size
    return
end

idx_prefix_size(kct::NeoKCT) = kct.idx[1].x

# WARNING: Will render "in place" KCT invalid (needs more work)
function collapse!(kct::NeoKCT{K, Ab, W, C}) where {K, Ab<:Alphabet, W<:Unsigned, C}
    deduped, perms, global_perms = permdedup(kct.counts)
    new_table = kct.table
    start_time = now()
    printstyled("Updating k-mer Indexes...", color=:green)
    @threads for (i, k_elem) in collect(enumerate(kct.table))
        new_table[i] = K_Element{K, Ab, C}(k_elem.seq, NTuple{length(k_elem.chunk_ids), UInt32}(getindex.(Ref(global_perms), k_elem.chunk_ids)))
    end
    printstyled("Table Collapse done in: $(now()-start_time)", color=:green)
    return NeoKCT{K, Ab, W, C}(new_table, deduped, kct.idx, kct.samples)
end

function Base.findfirst(kct::NeoKCT{K, Ab}, key::Mer{K, Ab}) where {K, Ab<:Alphabet}
    prefix_size = sizeof_idx_prefix(kct)
    symbol_size = bits_per_symbol(Ab())
    idx_key = key.data[1] >> (prefix_size*symbol_size) + 1 
    r = kct.idx[2][idx_key]
    i = searchsortedfirst(kct.table, key, r.start, r.stop, Base.Forward)
    return i != 0 && kct.table[i].seq == key ? i : 0
end

# Initializing a NeoKCT with the hashtable of the first sample
function NeoKCT{K, Ab, W}(sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet, W<:Unsigned}
    kct = NeoKCT{K, Ab, W}()
    @showprogress desc="parsing hash table into KCT..." for (k_bits, count) in sample_hashtable
        word_id = push!(kct.counts, UInt64(count))
        tmp_seq = Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,))
        push!(kct.table, K_Element{K, Ab}(tmp_seq, NTuple{1, UInt32}(UInt32(word_id))))
    end
    sort!(kct)
    compute_index!(kct)
    return kct
end

function _push_new_kmer_counts!(counts::PackedArray{UInt32, W}, prev_samples::Int, count::UInt32) where {W<:Unsigned}
    wid = new_word!(counts)  # allocate word, store nothing
    chunk_ids = UInt32[wid]
    @inbounds for _ in 1:prev_samples
        new_wid = push!(counts, W(0), wid)
        if new_wid != wid
            push!(chunk_ids, UInt32(new_wid))
            wid = new_wid
        end
    end
    new_wid = push!(counts, W(count & typemax(W)), wid)  # add new sample's count
    if new_wid != wid
        push!(chunk_ids, UInt32(new_wid))
    end
    return NTuple{length(chunk_ids), UInt32}(chunk_ids)
end

# Unecessary as trailing is assumed and accounted for in the push logic
function _pad_trailing_zeros!(kct::NeoKCT{K, Ab, W}, target_count::Int) where {K, Ab<:Alphabet, W<:Unsigned}
    @inbounds for i in eachindex(kct.table)
        ke = kct.table[i]
        total_stored = sum(sum(@view kct.counts.bitmap[word_bitmap_slice(Int(cid), W)]) for cid in ke.chunk_ids)
        last_wid = Int(ke.chunk_ids[end])
        for _ in 1:(target_count - total_stored)
            new_wid = push!(kct.counts, UInt64(0), last_wid)
            if new_wid != last_wid
                kct.table[i] = K_Element{K, Ab}(kct.table[i].seq, push(kct.table[i].chunk_ids, UInt32(new_wid)))
                last_wid = new_wid
            end
        end
    end
end

function find_shared_words(kct::NeoKCT)
    seen_wids = Set{UInt32}()
    shared_words = Set{UInt32}()
    @showprogress "Finding Shared Words..." for ke in kct.table
        for cid in ke.chunk_ids
            if cid in seen_wids
                push!(shared_words, cid)
            else
                push!(seen_wids, cid)
            end
        end
    end
    return shared_words
end


# Factimily push! function to interface NTuple like a vector for a more "painless" switch to them in chunk_ids
function push(tuple::NTuple{N, T}, val::T) where {N, T}
    return  NTuple{N+1, T}((tuple..., val))
end

# Main push! from a sample to an existing KCT
# TODO: This is getting very large, break up the logic across functions
function Base.push!(kct::NeoKCT{K, Ab, W}, sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet, W<:Unsigned}
    shared_words = find_shared_words(kct)

    @showprogress desc="Adding Sample $(kct.samples.x+1) to Table..." for (k_bits, count) in sample_hashtable
        tmp_seq = Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,))
        k_pos = findfirst(kct, tmp_seq)

        word_id = 0::Int
        if k_pos == 0
            new_chunk_ids = _push_new_kmer_counts!(kct.counts, kct.samples.x, count)
            ke = K_Element{K, Ab}(tmp_seq, new_chunk_ids)
            push!(kct.table, ke)
        else
            # Count how many values are already stored for this k-mer across all its chunk words.
            # If a k-mer was absent from intermediate samples, those zeros were never pushed,
            # so we must pad them now before appending the current sample's count.
            total_stored = sum(sum(@view kct.counts.bitmap[word_bitmap_slice(Int(cid), W)]) for cid in kct.table[k_pos].chunk_ids)
            missing_counts = kct.samples.x - total_stored
            last_wid = Int(kct.table[k_pos].chunk_ids[end])

            if last_wid in shared_words
                if missing_counts > 0
                    new_wid = push!(kct.counts, W(0)) # new word, first zero
                    kct.table[k_pos] = K_Element{K, Ab}(kct.table[k_pos].seq, push(kct.table[k_pos].chunk_ids, UInt32(new_wid)))
                    for _ in 2:missing_counts
                        next_wid = push!(kct.counts, W(0), new_wid)
                        if next_wid != new_wid
                            kct.table[k_pos] = K_Element{K, Ab}(kct.table[k_pos].seq, push(kct.table[k_pos].chunk_ids, UInt32(next_wid)))
                            new_wid = next_wid
                        end
                    end
                    word_id = push!(kct.counts, W(count & typemax(W)), new_wid)
                    if word_id != new_wid
                        kct.table[k_pos] = K_Element{K, Ab}(kct.table[k_pos].seq, push(kct.table[k_pos].chunk_ids, UInt32(word_id)))
                    end
                else
                    word_id = push!(kct.counts, W(count & typemax(W)))
                    kct.table[k_pos] = K_Element{K, Ab}(kct.table[k_pos].seq, push(kct.table[k_pos].chunk_ids, UInt32(word_id)))
                end
            else
                # Word is exclusively owned by this k-mer, safe to pack into it.
                for _ in 1:missing_counts
                    new_wid = push!(kct.counts, W(0), last_wid)
                    if last_wid != new_wid
                        kct.table[k_pos] = K_Element{K, Ab}(kct.table[k_pos].seq, push(kct.table[k_pos].chunk_ids, UInt32(new_wid)))
                        last_wid = new_wid
                    end
                end

                word_id = push!(kct.counts, W(count & typemax(W)), last_wid)
                if last_wid != word_id
                    kct.table[k_pos] = K_Element{K, Ab}(kct.table[k_pos].seq, push(kct.table[k_pos].chunk_ids, UInt32(word_id)))
                end
            end
        end
    end

    sort!(kct)
    compute_index!(kct)
    # _pad_trailing_zeros!(kct, kct.samples.x + 1)  # pad to new sample count
    kct.samples.x += 1
    # return kct
end

function assemble_count_vector(kct::NeoKCT{K, Ab, W}, chunk_ids::NTuple) where {K, Ab<:Alphabet, W<:Unsigned}
    counts = reduce(vcat, [assemble_word(kct.counts, i) for i in chunk_ids])
    n_missing = kct.samples.x - length(counts)
    return n_missing > 0 ? vcat(counts, zeros(W, n_missing)) : counts
end

function Base.getindex(kct::NeoKCT{K, Ab, W}, i::Integer)  where {K, Ab<:Alphabet, W<:Unsigned}
    return kct.table[i].seq => assemble_count_vector(kct, kct.table[i].chunk_ids)
end

function build_kct(sample::String, K::Int=30, chunks::Int = 500_000; word_size::DataType=UInt64, collapse::Bool=true)
    kct = NeoKCT{K÷3, AAAlphabet, word_size}(jello_superthreaded_hash(sample, K, chunks))
    kct = collapse ? collapse!(kct) : kct
    return kct
end

function build_kct(samples::AbstractVector{String}, K::Int=30, chunks::Int = 500_000; word_size::DataType=UInt64, save_at_samples::AbstractVector{Int}=Int[], save_path::String = "", collapse_every::Int=1, benchmark_every::Int=5)
    kct = NeoKCT{K÷3, AAAlphabet, word_size}(jello_superthreaded_hash(popfirst!(samples), K, chunks))
    mkpath(save_path*"benchmarks/")
    for (i, sample) in enumerate(samples)
        sample_hashtable = jello_superthreaded_hash(sample, K, chunks)
        push!(kct, sample_hashtable)
        kct = ((i+1) % collapse_every == 0) ? collapse!(kct) : kct 
        i+1 in save_at_samples && write_kct(kct, save_path*"$(day(now()))_$(month(now()))_$(year(now()))_[$(i+1)_samples]_neokct.kct")
        ((i+1) % benchmark_every == 0) && benchmark_kct(kct, save_path*"benchmarks/")
    end
    return kct
end

### EXPERIMENTAL

# temporary path to quickly test on all samples
samples = readlines(open("/u/jacquinn/phd_stuff/data/tcga_fastqs.paths", "r"))

function benchmark_kct(kct::NeoKCT{K, Ab}, benchmark_path::String) where {K, Ab<:Alphabet}

    # Benchmark kct element sizes 
    f = Figure()

    xs = ["K-mers", "Indexes", "Counts", "Bitmap"]
    ys = [
        Base.summarysize(kct.table[begin].seq) * length(kct.table),  # K-mer sequences are constant sizes
        sum([Base.summarysize(k_elem.chunk_ids) for k_elem in kct.table]),  # Indexes are variable and need to be measured systematically
        Base.summarysize(kct.counts.words),
        Base.summarysize(kct.counts.bitmap)
    ]

    Axis(f[1, 1],
         title = "$(kct.samples.x) Samples NeoKCT - Component Sizes",
         subtitle = "Total Size: $(Base.format_bytes(sum(ys))) - [Prefix Index Cost: $(Base.format_bytes(Base.summarysize(kct.idx)))]",
         ylabel = "Size (Bytes)",
         xticks = (1:length(xs), xs)
         )

    barplot!(ys, 
             color = ys, strokecolor = :black, strokewidth = 1,
             bar_labels = ys)

    CairoMakie.save(benchmark_path*"sizes_benchmark_$(kct.samples.x)_samples.svg", f)

    # Benchmark kct k-mer request speed
    k_mers = getfield.(rand(kct.table, 10_000_000), :seq)
    start = now()
    for k_mer in k_mers
        findfirst(k->k==k_mer, kct)
    end
    query_time = now()-start

    # Benchmark count table


    # Benchmark count indexes

    # Read Benchmark file

    # Compute elapsed time and plot

    # Update benchmark file

end

function p_sort_merge(kct::NeoKCT{K, Ab, W}, sample_hashtable::Dict{UInt64, UInt32}, n_chunks::Int=200) where {K, Ab<:Alphabet, W<:Unsigned}
    right = psort!(collect(keys(sample_hashtable))) # bit format
    left = kct.table  # K_Element format (use .seq.data[1] to access bits)
    
    cutoffs_left = length(left)÷n_chunks:length(left)÷n_chunks:length(left)
    cutoffs_right = searchsortedfirst.(Ref(right), left[cutoffs_left])

    merge_blocks = Vector{Vector{K_Element{K, Ab}}}(undef, n_chunks)
    threads = Vector{Task}(undef, n_chunks)
    prev_left = 0
    prev_right = 0

    @showprogress for thread in 1:nchunks
        left_block = left[prev_left:cutoffs_left[thread]]
        right_block = right[prev_right:cutoffs_right[thread]]
        threads[thread] = @spawn merge_blocks[i] = smerge(left_bloc, right_block)

        prev_left, prev_right = cutoffs_left[thread], cutoffs_right[thread]
    end

    wait.(threads)
    return ## TODO: Combine sub K_Element table
end


include("KCTLoader.jl")