if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(".")
end

using Kmers
using BioSequences
using BioSymbols
using ProgressMeter
using Base.Threads

include("BioParser.jl")
include("AAAlphabet.jl")

# Parallel k-mer counter. Streams a FASTQ (or other bio file) in chunks, counts each chunk
# on its own thread as amino-acid k-mers, and folds the per-chunk hash-tables together as
# they finish. The entry point is jello_superthreaded_hash.

# 2-bit code (0..3) for a nucleic acid.
function twobit(n::NucleicAcid)::UInt8
    bits = compatbits(n)
    return UInt8(trailing_zeros(bits))
end

# Translate a DNA/RNA k-mer in a single reading frame to a 5-bit amino-acid k-mer.
# Returns nothing if an in-frame stop codon is hit or the length is not a multiple of 3.
function translate(kmer::Kmer{Ab}, code::BioSequences.GeneticCode=ncbi_trans_table[1]) where {Ab<:Union{DNAAlphabet, RNAAlphabet}}
    bps = BioSequences.bits_per_symbol(BioSequences.BitsPerSymbol(AAAlphabet()))
    aa_sequence = AminoAcid[]
    for interval in Iterators.partition(1:length(kmer), 3)
        length(interval) != 3 && break
        code_idx = UInt64(0)
        for (i, k_i) in enumerate(interval)
            code_idx |= twobit(kmer[k_i]) << (4-(i-1)*2)
        end
        aa = code[code_idx]
        aa == AA_Term && return
        push!(aa_sequence, aa)
        
    end
    KP = length(aa_sequence)
    CP = Int64(ceil(bps*KP/64))
    pkmer = Kmer{AAAlphabet, KP, CP}(aa_sequence)
    return pkmer
end

# Chop a DNA/RNA sequence into its overlapping K-mers.
function k_merize(sequence::LongSequence{Ab}; K::Int) where {Ab<:Alphabet}
    length(sequence) < K && error("Read of size $(length(sequence)) too small for $K-mers")
    to_return = Kmer{Ab, K}[]
    for i in 1:(length(sequence)-K+1)
        push!(to_return, Kmer{Ab, K}(sequence[i:i+K-1]))
    end
    return to_return
end

# Chop a 5-bit amino-acid sequence into its overlapping K-mers.
function k_merize(sequence::LongAA; K::Int)
    length(sequence) < K && error("Read of size $(length(sequence)) too small for $K-mers")
    to_return = Kmer{AAAlphabet, K}[]
    for i in 1:(length(sequence)-K+1)
        push!(to_return, Kmer{AAAlphabet, K}(sequence[i:i+K-1]))
    end
    return to_return
end

# Read through a bio file (see BioParser), filling fixed-size chunks of sequences. Each full
# chunk is handed to a spawned task for k-mer counting, and the leftover partial chunk is
# flushed at EOF. Blocks until every counting task has finished.
function chunk_stream(file::String, K::Int, merge_queue::Channel{Dict{UInt64, UInt32}}, chunking::Int=1_000_000; verbose::Bool=false)
    counting_tasks = Vector{Task}()
    current_chunk = String[]
    sizehint!(current_chunk, chunking)
    i = 0
    io = stream(file)

    progress = ProgressUnknown(desc = "Processing $chunking k-mers chunks...")
    for seq in io
        push!(current_chunk, seq)
        i += 1
        if i % chunking == 0
            let chunk = current_chunk
                push!(counting_tasks, @spawn count_kmers(chunk, K, merge_queue, verbose=true))
            end
            current_chunk = String[]
            sizehint!(current_chunk, chunking)

        end

        !verbose && next!(progress; showvalues=[
            ("active counting tasks", length(filter(x->x.state!=:done, counting_tasks))),
            ("items in merge queue", merge_queue.n_avail_items),
            ("chunks sent", i÷chunking)
        ])
    end

    if !isempty(current_chunk)
        push!(counting_tasks, @spawn count_kmers(current_chunk, K, merge_queue, verbose=true))
        !verbose && ProgressMeter.update!(progress; showvalues=[
            ("active counting tasks", length(filter(x->x.state!=:done, counting_tasks))),
            ("items in merge queue", merge_queue.n_avail_items),
            ("chunks sent", i÷chunking+1)
        ])
    end
    wait.(counting_tasks)
    finish!(progress)
end

# Count one chunk of DNA/RNA sequences (from chunk_stream) as translated 5-bit AA k-mers,
# then drop the resulting hash-table onto merge_queue. Reads with an 'N' are skipped. Any
# error is logged and the (possibly partial) hash-table is still enqueued.
function count_kmers(chunk::Vector{String}, K::Int, merge_queue::Channel{Dict{UInt64, UInt32}}; verbose::Bool=false)
    hash = Dict{UInt64, UInt32}()
    try 
        for l in chunk
            'N' in l && continue
            seq = LongSequence{DNAAlphabet{2}}(l)
            length(seq) < K && continue
            kmers = translate.(k_merize(seq, K=K))
            for kmer in kmers
                isnothing(kmer) && continue
                if !(haskey(hash, kmer.data[1]))
                    hash[kmer.data[1]] = UInt32(0)
                end
                hash[kmer.data[1]] += 1
            end
        end
    catch e
        @error "counting thread task crashed during hash building" exception=(e, catch_backtrace())
    end
    try
        put!(merge_queue, hash)
    catch e
        @error "counting thread task crashed during merge queue put" exception=(e, catch_backtrace())
    end
end

"""
    jello_superthreaded_hash(fastq, K, chunking=1_000_000, queue_size=128; verbose=false)
        -> Dict{UInt64, UInt32}

Count every amino-acid `K÷3`-mer in `fastq` and return a map from a k-mer's raw
bit-encoding (`Kmer{AAAlphabet, K÷3, 1}.data[1]`) to its count. `K` is the
nucleotide k-mer length. Reads are k-merized, translated in a single frame, and
counted per chunk on separate threads, while a pairer task keeps merging the
finished per-chunk hash-tables in a binary-tree pattern until one table is left.
`chunking` sets the reads per chunk, `queue_size` the merge-queue depth. This is
the per-sample input `KCT{K÷3, AAAlphabet}` and `push!` expect.
"""
function jello_superthreaded_hash(fastq::String, K::Int, chunking::Int=1_000_000, queue_size::Int=128; verbose::Bool=false)
    Hash_Type = Dict{UInt64, UInt32}
    merge_queue = Channel{Hash_Type}(queue_size)
    merge_tasks = Task[]
    paired = Hash_Type[]
    signal = Channel{Nothing}(Inf)
    inflight = Threads.Atomic{Int64}(0)

    pairer = @async try
        for child_hash in merge_queue
            push!(paired, child_hash)
            if length(paired) >= 2
                left, right = pop!(paired), pop!(paired)
                task = @spawn begin
                    Threads.atomic_add!(inflight, 1)
                    put!(merge_queue, merge(+, left, right))
                    Threads.atomic_sub!(inflight, 1)
                    put!(signal, nothing)
                end
                push!(merge_tasks, task)
            end
        end
    catch e
        @error "Pairer task crashed" exception=(e, catch_backtrace())
    end

    chunk_stream(fastq, K, merge_queue, chunking, verbose = verbose)

    progress = ProgressUnknown(desc = "Counting tasks done. Waiting on merger tasks to complete...")
    while inflight[] != 0 || merge_queue.n_avail_items != 0 || length(paired) != 1
        take!(signal)
        !verbose && next!(progress; showvalues=[
            ("active merging tasks", length(filter(x->x.state!=:done, merge_tasks))),
            ("items in merge queue", merge_queue.n_avail_items),
        ])
    end
    close(merge_queue)
    wait(pairer)
    finish!(progress)
    return pop!(paired)
end
