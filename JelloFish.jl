using BioSequences
using BioSymbols
using Kmers
using ProgressMeter
using Base.Threads
using EzXML
using GZip

include("BioParser.jl")
include("AAAlphabet.jl")

# Returns 2 bit version of a nucleic acid
function twobit(n::NucleicAcid)::UInt8
    bits = compatbits(n)
    return UInt8(trailing_zeros(bits))
end

# Translates a K-mer of DNA/RNA to a 5 bits amino acid k-mer
function translate(kmer::Kmer{Ab}, code::BioSequences.GeneticCode=ncbi_trans_table[1]) where {Ab<:Union{DNAAlphabet, RNAAlphabet}}
    bps = BioSequences.bits_per_symbol(BioSequences.BitsPerSymbol(AAAlphabet()))
    aa_sequence = AminoAcid[]
    for interval in Iterators.partition(1:length(kmer), 3)
        length(interval) != 3 && break
        code_idx = UInt64(0)
        for (i, k_i) in enumerate(interval)
            # kmer[k_i] == DNA_N && return
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

# Chops a DNA/RNA sequence into an array of K-mers 
function k_merize(sequence::LongSequence{Ab}; K::Int) where {Ab<:Alphabet}
    length(sequence) < K && error("Read of size $(length(sequence)) too small for $K-mers")
    to_return = Kmer{Ab, K}[]
    for i in 1:(length(sequence)-K+1)
        push!(to_return, Kmer{Ab, K}(sequence[i:i+K-1]))
    end
    return to_return
end

# Chops a 5 bits amino acid sequence into an array of K-mers 
function k_merize(sequence::LongAA; K::Int)
    length(sequence) < K && error("Read of size $(length(sequence)) too small for $K-mers")
    to_return = Kmer{AAAlphabet, K}[]
    for i in 1:(length(sequence)-K+1)
        push!(to_return, Kmer{AAAlphabet, K}(sequence[i:i+K-1]))
    end
    return to_return
end

# Reads through a bio file (see BioParser), populating chunks of Sequences
# Full chunks are dispatched to a thread for k-mer counting
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
            to_process = current_chunk
            push!(counting_tasks, @spawn count_kmers(to_process, K, merge_queue,  verbose=true))
            current_chunk = String[]
            sizehint!(current_chunk, chunking)

        end
        !verbose && next!(progress; showvalues=[
            ("active counting tasks", length(filter(x->x.state!=:done, counting_tasks))),
            ("items in merge queue", merge_queue.n_avail_items),
            ("chunks sent", i÷chunking)
        ])
        # i >= 5_000_000 && break  ## REMOVE ME WHEN DONE TESTING
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

# K-mer counting (with translation to 5 bits AA) of a chunk of DNA/RNA sequences (from chunk_stream)
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

# Main function of JelloFish. Spawns chunk_stream task to read a fastq, continuously merges count hash-tables until done
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
                    start = now()
                    Threads.atomic_add!(inflight, 1)
                    put!(merge_queue, merge(+, left, right))
                    Threads.atomic_sub!(inflight, 1)
                    # verbose && tprintln("Completed merging in $(now()-start)")
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
