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
# include("JaggedRLEArrays.jl")
include("PackedArray.jl")
include("AAAlphabet.jl")
include("BidirArray.jl")
include("JelloFish.jl")

global const VERSION = 1.0

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

Base.isless(a::K_Element, b::K_Element) = a.seq < b.seq 
Base.isless(a::Kmer, b::K_Element) = a < b.seq
Base.isless(a::K_Element, b::Kmer) = a.seq < b
Base.:(==)(a::K_Element, b::K_Element) = a.seq == b.seq


"""
A sorted k-mer count table implementing bitpacked count arrays.

    NeoKCT{K, Ab<:Alphabet, W<:Unsigned}()
"""
struct NeoKCT{K, Ab<:Alphabet, W<:Unsigned, C}
    table::Vector{K_Element{K, Ab, C}}
    counts::PackedArray{UInt32, W}
    idx::Vector{UnitRange{Int64}}
    samples::Base.RefValue{Int64}  # Send me to jail for this one
    version::Float64
end

NeoKCT{K, Ab, W, C}(table, counts, idx, samples) where {K, Ab<:Alphabet, W<:Unsigned, C} = NeoKCT(table, counts, idx, samples, VERSION)
NeoKCT{K, Ab, W}() where {K, Ab<:Alphabet, W<:Unsigned} = NeoKCT(K_Element{K, Ab, 1}[], PackedArray{UInt32, W}(), fill(0:-1, 4^15), Ref(1), VERSION)

Base.sort!(kct::NeoKCT) = psort!(kct.table)

function compute_index!(kct::NeoKCT{K, Ab}) where {K, Ab<:Alphabet}
    start = 1
    last_key = 0x0000000000000000
    symbol_size = bits_per_symbol(Ab())
    # println(symbol_size)
    
    @showprogress "Computing Binary Search Index..." for i in eachindex(kct.table)
        key = kct.table[i].seq.data[1] >> (Int(round(K*0.4))*symbol_size)
        if key > last_key
            # println(last_key)
            # last_key+1 > length(kct.idx) && append!(kct.idx, fill(0:-1, length(kct.idx)))  # Currently memory leaks
            kct.idx[last_key+1] = start:i
            start = i
            last_key = key
        end
    end
    kct.idx[last_key+1] = start:length(kct.table)
    return
end

#IDEA: Generated function collapse

# WARNING: Will render "in place" KCT invalid (needs more work)
function collapse!(kct::NeoKCT{K, Ab, W, C}) where {K, Ab<:Alphabet, W<:Unsigned, C}
    deduped, perms, global_perms = permdedup(kct.counts)
    new_table = kct.table
    start_time = now()
    printstyled("Updating k-mer Indexes...", color=:green)
    @threads for (i, k_elem) in collect(enumerate(kct.table))
        new_table[i] = K_Element{K, Ab, C}(k_elem.seq, NTuple{length(k_elem.chunk_ids), UInt32}(getindex.(Ref(global_perms), k_elem.chunk_ids)))
        # for (j, prev_idx) in enumerate(k_elem.chunk_ids)
        #     new_table[i].chunk_ids[j] = global_perms[prev_idx]
        # end
    end
    printstyled("Table Collapse done in: $(now()-start_time)", color=:green)
    return NeoKCT{K, Ab, W, C}(new_table, deduped, kct.idx, kct.samples)
end

# function Base.findfirst(kct::NeoKCT{K, Ab}, key::Mer{K, Ab}) where {K, Ab<:Alphabet}
#     symbol_size = bits_per_symbol(Ab())
#     idx_key = key.data[1] >> (Int(round(K*0.35))*symbol_size) + 1 
#     r = kct.idx[idx_key]
#     # r = 1:length(kct.table)
#     t = @view kct.table[r]
#     i = searchsortedfirst(t, key)
#     if i > length(t) || t[i].seq != key
#         return 0
#     end
#     return r[1] + i-1
# end

function Base.findfirst(kct::NeoKCT{K, Ab}, key::Mer{K, Ab}) where {K, Ab<:Alphabet}
    symbol_size = bits_per_symbol(Ab())
    idx_key = key.data[1] >> (Int(round(K*0.4))*symbol_size) + 1 
    r = kct.idx[idx_key]
    # println(r[1], r[2])
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
    # kct = collapse ? collapse!(kct) : kct
    return kct
end

# Pushing a new sample to the NeoKCT
# function Base.push!(kct::NeoKCT{K, Ab}, sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet}
#     @showprogress desc="Adding sample $(kct.samples+1) to table" for (k_bits, count) in sample_hashtable
#         tmp_seq = Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,))
#         k_pos = findfi0rst(kct, tmp_seq)

#         word_id = k_pos == 0 ?
#             push!.(Ref(kct.counts), UInt64.(vcat(zeros(kct.samples), count))) :
#             push!(kct.counts, UInt64(count), kct.table[k_pos].chunk_ids[end])
#         k_pos = k_pos == 0 ? length(push!(kct.table, K_Element{K, Ab}(tmp_seq, [word_id]))) : k_pos
#         if kct.table[k_pos].chunk_ids[end] != word_id
#             push!(kct.table[k_pos].chunk_ids, word_id)
#         end
#     end
#     sort!(kct)
#     compute_index!(kct)
#     return kct
# end

# function Base.push!(kct::NeoKCT{K, Ab}, sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet}
#     @showprogress desc="Adding sample $(kct.samples+1) to table" for (k_bits, count) in sample_hashtable
#         tmp_seq = Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,))
#         k_pos = findfirst(kct, tmp_seq)

#         word_id = k_pos == 0 ?
#             push!.(Ref(kct.counts), UInt64.(vcat(zeros(kct.samples), count))) :
#             push!(kct.counts, UInt64(count), k_pos)
#         k_pos = k_pos == 0 ? length(push!(kct.table, K_Element{K, Ab}(tmp_seq, [word_id]))) : k_pos
#         push!(kct.table[k_pos].chunk_ids, word_id)
#     end
#     sort!(kct)
#     compute_index!(kct)
#     return kct
# end

function _push_new_kmer_counts!(counts::PackedArray{UInt32,W}, prev_samples::Int, count::UInt32) where {W<:Unsigned}
    # Start a new packed word by pushing one zero
    wid = push!(counts, UInt64(0))

    # Pad zeros for the previous samples
    @inbounds for _ in 2:prev_samples
        wid = push!(counts, UInt64(0), wid)
    end

    # Push the count for this new sample
    wid = push!(counts, UInt64(count), wid)
    return wid
end

function Base.push!(kct::NeoKCT{K, Ab, W}, sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet, W<:Unsigned}
    @showprogress desc="Adding Sample $(kct.samples.x+1) to Table..." for (k_bits, count) in sample_hashtable
        tmp_seq = Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,))
        k_pos = findfirst(kct, tmp_seq)

        word_id = 0::Int
        if k_pos == 0
            word_id = _push_new_kmer_counts!(kct.counts, kct.samples.x, count)

            ke = K_Element{K, Ab}(tmp_seq, NTuple{1, UInt32}(UInt32(word_id)))
            push!(kct.table, ke)
        else
            word_id = push!(kct.counts, UInt64(count), kct.table[k_pos].chunk_ids[end])

            if kct.table[k_pos].chunk_ids[end] != UInt32(word_id)
                kct.table[k_pos] = K_Element{K, Ab}(kct.table[k_pos].seq, push(kct.table[k_pos].chunk_ids, UInt32(word_id)))
            end
        end
    end

    sort!(kct)
    compute_index!(kct)
    
    kct.samples.x += 1
    # return kct
end

function assemble_count_vector(kct::NeoKCT{K, Ab, W}, chunk_ids::NTuple) where {K, Ab<:Alphabet, W<:Unsigned}
    return reduce(vcat, [assemble_word(kct.counts, i) for i in chunk_ids])
end

function Base.getindex(kct::NeoKCT{K, Ab, W}, i::Integer)  where {K, Ab<:Alphabet, W<:Unsigned}
    return kct.table[i].seq => assemble_count_vector(kct, kct.table[i].chunk_ids)
end

### EXPERIMENTAL

# temporary path to quickly test on all samples
samples = readlines(open("/u/jacquinn/phd_stuff/data/tcga_fastqs.paths", "r"))


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
        i+1 in save_at_samples && save(kct, save_path*"$(day(now()))_$(month(now()))_$(year(now()))_[$(i+1)_samples]_neokct.kct")
        ((i+1) % benchmark_every == 0) && benchmark_kct(kct, save_path*"benchmarks/")
    end
    return kct
end

to_type(s::String) = eval(Symbol(s))

function save(kct::NeoKCT{K, Ab}, fn::String) where {K, Ab<:Alphabet}
    f = open(fn, "w")
    write(fn, kct)
    close(f)
end

function Base.write(io::IO, tuple::Tuple)
    write(io, length(tuple))
    write.(io, tuple)
end

function Base.write(io::IO, k_elem::K_Element{K, Ab, C}) where {K, Ab<:AAAlphabet, C}
    write(io, K)
    write(io, codeunits(String(Ab.name.singletonname)))
    write(io, C)
    write.(io, k_elem.seq.data)
    write(io, k_elem.chunk_ids)
end

function Base.write(io::IO, kct::NeoKCT{K, Ab}) where {K, Ab<:Alphabet}
    write(io, kct.version)
    write(io, K)
    write(io, (String(Ab.name.singletonname)))
    write(io, length(kct.table))
    write(io, length(kct.counts.words))
    write(io, length(kct.counts.bitmap))
    write.(io, kct.table)
    write(io, kct.counts.words)
    write(io, kct.counts.bitmap)
    write(io, kct.samples.x)
end

function load(path::String)
    io = open(path, "r")
    version = read(io, Float64)
    K = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:10])))  # TODO: Fix this abomination
    table_length = read(io, Int64)
    words_length = read(io, Int64)
    bitmap_length = read(io, Int64)
    table = Vector{K_Element{K, Ab}}(undef, table_length)
    @showprogress for i in 1:table_length
        K = read(io, Int64)
        Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:10])))  # TODO: Fix this abomination
        kmer_C = read(io, Int64)
        data = Tuple(read(io, UInt64) for _ in 1:kmer_C)
        word_C = read(io, Int64)
        chunk_ids = Tuple(read(io, UInt32) for _ in 1:word_C)
        table[i] = K_Element{K, Ab, kmer_C}(Kmer{Ab, K, kmer_C}(Kmers.unsafe, data), chunk_ids)
    end
    words = Vector{UInt64}(undef, words_length)
    @showprogress for i in 1:words_length
        words[i] = read(io, UInt64)
    end
    bitmap = BitVector(falses(bitmap_length))
    @showprogress for i in 1:bitmap_length  ## DIVIDE BY WORD SIZE, READ BY WORD SIZE
        bitmap[i] = read(io, Bool)
    end
    samples = Ref(read(io, Int64))


end

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
         subtitle = "Total Size: $(Base.format_bytes(sum(ys)))",
         ylabel = "Size (Bytes)",
         xticks = (1:length(xs), xs)
         )

    barplot!(ys, 
             color = ys, strokecolor = :black, strokewidth = 1,
             bar_labels = ys)

    CairoMakie.save(benchmark_path*"sizes_benchmark_$(kct.samples.x)_samples.svg", f)

    # Benchmark kct k-mer request speed
    # k_mers = getfield.(rand(kct.table, 10_000_000), :seq)
    # start = now()
    # for k_mer in k_mers
    #     findfirst(k->k==k_mer, kct)
    # end
    # query_time = now()-start


    # Benchmark count table

    # Benchmark count indexes

end

# function sort_p_merge

# get_bits(a::UInt) = a
# get_bits(a::K_Element) = a.seq.data[1]

# Base.isless(a::Union{K_Element, UInt}, b::Union{K_Element, UInt}) = get_bits(a) < get_bits(b) 
# Base.:(==)(a::Union{K_Element, UInt}, b::Union{K_Element, UInt}) = get_bits(a) == get_bits(b) 

Base.isless(a::K_Element, b::UInt) = a.seq.data[1] < b 
Base.isless(a::UInt, b::K_Element) = a < b.seq.data[1]
Base.:(==)(a::K_Element, b::UInt) = a.seq.data[1] == b
Base.:(==)(a::UInt, b::K_Element) = a == b.seq.data[1]


# Factimily push! function to interface NTuple like a vector for a more "painless" switch to them in chunk_ids
function push(tuple::NTuple{N, T}, val::T) where {N, T}
    return  NTuple{N+1, T}((tuple..., val))
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

function testing()
    times = []
    fastq = "/u/jacquinn/Aboleth/data/gdc_fastq_tests/COAD_TCGA-AA-3688-01A-01R-0905-07/100908_UNC7-RDR3001641_00027_FC_62ER1AAXX.6_R1.fastq"
    start = now()
    result_sc = jello_superthreaded_hash(fastq, 30, 500_000)
    kct = NeoKCT{10, AAAlphabet, UInt64}(result_sc)
    kct = collapse(kct)
    push!(times, now()-start)
    for i in 1:16
        start = now()
        result_sc = jello_superthreaded_hash(fastq, 30, 500_000)
        push!(kct, result_sc)
        kct = collapse(kct)
        push!(times, now()-start)
    end
    println(times)
end

function dedup_from_perm(vals::AbstractVector, perm::AbstractVector{<:Integer})
    n = length(perm)
    unique_vals = Vector{eltype(vals)}()
    unique_perm = Int[]
    invmap = similar(perm, Int)

    if n == 0
        return unique_vals, unique_perm, invmap
    end

    # Work on the sorted order
    first_idx = perm[1]
    last_val = vals[first_idx]
    push!(unique_vals, last_val)
    push!(unique_perm, first_idx)
    current_unique = 1
    invmap[first_idx] = current_unique

    @inbounds for k in 2:n
        idx = perm[k]
        v = vals[idx]
        if v != last_val
            last_val = v
            current_unique += 1
            push!(unique_vals, v)
            push!(unique_perm, idx)
        end
        invmap[idx] = current_unique
    end

    return unique_vals, unique_perm, invmap
end