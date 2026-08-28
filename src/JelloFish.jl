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

## Fast rolling counter ##

#= The `translate` + `k_merize` path above allocates a Vector{Kmer} per read and constructs
a fresh Kmer per position, which at ~1e10 read-k-mers per sample is the whole runtime. The
counter below rolls a 2-bit code one base at a time with nothing allocated in the inner
loop. `translate=false` emits the raw nucleotide K-mer code, `translate=true` folds the
K-nt window through a 64-entry codon table into a 5*(K/3)-bit AA code. AA output is
bit-identical to `translate(k_merize(...))` (same codon table, same stop-drop, same symbol
order), so tables built either way stay comparable. =#

# ASCII byte -> 2-bit code (A/a=0 C/c=1 G/g=2 T/t=3), 0xff for anything else.
const _NT2BIT = let t = fill(0xff, 256)
    for (c, v) in (('A', 0x00), ('C', 0x01), ('G', 0x02), ('T', 0x03))
        t[UInt8(c) + 1] = v
        t[UInt8(lowercase(c)) + 1] = v
    end
    t
end

# 6-bit codon (b1<<4 | b2<<2 | b3, same 2-bit order as _NT2BIT) -> 5-bit AAAlphabet code,
# 0xff on a stop codon. Built from the exact GeneticCode indexing the `translate` above uses.
const _CODON2AA5 = let t = Vector{UInt8}(undef, 64), gc = BioSequences.ncbi_trans_table[1]
    for c in 0:63
        aa = gc[UInt64(c)]
        t[c + 1] = aa == AA_Term ? 0xff : UInt8(BioSequences.encode(AAAlphabet(), aa))
    end
    t
end

const _STOP_KEY = typemax(UInt64)  # sentinel from _fold_codons for an in-frame stop

# Fold the K-nt window packed in the low 2K bits of `code` (5' base highest) into a
# 5*kk-bit AA code, first codon in the high bits. _STOP_KEY if any codon is a stop.
@inline function _fold_codons(code::UInt64, twoK::Int, kk::Int)
    aac = UInt64(0)
    @inbounds for j in 0:(kk - 1)
        cod = Int((code >> (twoK - 6 - 6j)) & 0x3f)
        a = _CODON2AA5[cod + 1]
        a == 0xff && return _STOP_KEY
        aac = (aac << 5) | a
    end
    return aac
end

# True if the read has any byte that is not A/C/G/T. Matches the old "skip read with N",
# and also catches other IUPAC codes that used to crash the whole chunk in LongSequence().
@inline function _has_nonacgt(l::String)
    @inbounds for p in 1:ncodeunits(l)
        _NT2BIT[codeunit(l, p) + 1] == 0xff && return true
    end
    return false
end

# Read through a bio file (see BioParser), filling fixed-size chunks of reads. Each full
# chunk is handed to a spawned task for counting; the leftover partial chunk is flushed at
# EOF. `max_inflight` caps concurrently-running count tasks: the reader blocks on the
# semaphore once that many are outstanding, so ingestion cannot outrun counting. Blocks
# until every counting task has finished.
function chunk_stream(file::String, K::Int, merge_queue::Channel{Dict{UInt64, UInt32}},
                      chunking::Int=1_000_000; translate::Bool=false,
                      max_inflight::Int=max(2, Threads.nthreads()), verbose::Bool=false)
    counting_tasks = Task[]
    gate = Base.Semaphore(max_inflight)
    inflight = Threads.Atomic{Int}(0)
    current_chunk = String[]
    sizehint!(current_chunk, chunking)
    i = 0
    io = stream(file)

    spawn_chunk! = function (chunk)
        Base.acquire(gate)
        Threads.atomic_add!(inflight, 1)
        push!(counting_tasks, Threads.@spawn begin
            try
                count_kmers(chunk, K, merge_queue; translate = translate)
            finally
                Threads.atomic_sub!(inflight, 1)
                Base.release(gate)
            end
        end)
    end

    progress = ProgressUnknown(desc = "Processing $chunking read chunks...")
    for seq in io
        push!(current_chunk, seq)
        i += 1
        if i % chunking == 0
            let chunk = current_chunk
                spawn_chunk!(chunk)
            end
            current_chunk = String[]
            sizehint!(current_chunk, chunking)
        end
        !verbose && next!(progress; showvalues = [
            ("running counting tasks", inflight[]),
            ("items in merge queue", merge_queue.n_avail_items),
            ("chunks sent", i ÷ chunking),
        ])
    end

    isempty(current_chunk) || spawn_chunk!(current_chunk)
    wait.(counting_tasks)
    finish!(progress)
end

# Count one chunk of reads and drop the hash-table onto merge_queue. `translate=false`
# counts raw nucleotide K-mers, `translate=true` counts their in-frame AA (K/3)-mers.
# Reads with any non-ACGT base are skipped whole. Errors are logged and the partial table
# is still enqueued.
function count_kmers(chunk::Vector{String}, K::Int, merge_queue::Channel{Dict{UInt64, UInt32}};
                     translate::Bool=false, verbose::Bool=false)
    hash = Dict{UInt64, UInt32}()
    twoK = 2K
    twoK <= 64 || error("K=$K too large for a 64-bit nucleotide code")
    translate && (5 * (K ÷ 3) <= 64 || error("K=$K too large for a 64-bit AA code"))
    kmer_mask = twoK == 64 ? typemax(UInt64) : (UInt64(1) << twoK) - UInt64(1)
    kk = translate ? K ÷ 3 : K
    try
        @inbounds for l in chunk
            (ncodeunits(l) < K || _has_nonacgt(l)) && continue
            code = UInt64(0)
            valid = 0
            for p in 1:ncodeunits(l)
                code = ((code << 2) | _NT2BIT[codeunit(l, p) + 1]) & kmer_mask
                valid += 1
                valid < K && continue
                if translate
                    key = _fold_codons(code, twoK, kk)
                    key == _STOP_KEY && continue
                else
                    key = code
                end
                hash[key] = get(hash, key, UInt32(0)) + UInt32(1)
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
    jello_superthreaded_hash(fastq, K, chunking=1_000_000, queue_size=128;
                             translate=false, max_inflight=nthreads(), verbose=false)
        -> Dict{UInt64, UInt32}

Count k-mers in `fastq` and return a map from a k-mer's raw bit-encoding to its count.
`K` is always the nucleotide k-mer length. With `translate=false` (default) the keys are
raw `K`-nt codes (`Kmer{DNAAlphabet{2}, K, 1}.data[1]`); with `translate=true` they are the
in-frame amino-acid `K÷3`-mer codes (`Kmer{AAAlphabet, K÷3, 1}.data[1]`), stop codons
dropped. Reads are counted per chunk on separate threads while a pairer task folds the
finished per-chunk tables together in a binary tree. `chunking` sets reads per chunk,
`queue_size` the merge-queue depth, `max_inflight` the cap on concurrent counting tasks so
the reader cannot outrun the counters.
"""
function jello_superthreaded_hash(fastq::String, K::Int, chunking::Int=1_000_000, queue_size::Int=128;
                                  translate::Bool=false, max_inflight::Int=max(2, Threads.nthreads()),
                                  verbose::Bool=false)
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

    chunk_stream(fastq, K, merge_queue, chunking; translate = translate,
                 max_inflight = max_inflight, verbose = verbose)

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
