using BioSequences
using BioSymbols
using Kmers
using ProgressMeter
using Base.Threads
using EzXML
using GZip

# include("JaggedRLEArrays.jl")
# include("RLEKct.jl")
# include("thread_utils.jl")

struct Node{Sb<:BioSymbol}
    value::Sb
    nexts::Vector{Node{Sb}}
    count_idx::UInt32
end

struct Trie{Sb<:BioSymbol}
    K::Int
    roots::Vector{Node{Sb}}
    counts::Vector{UInt32}
end

function Base.push!(trie::Trie{Sb}, kmer::Km) where {Sb<:BioSymbol, Km<:Kmer}
    nexts = trie.roots
    last = nothing
    for (i, n) in enumerate(kmer)
        next = get_next(nexts, n)
        # println(next)
        if isnothing(next)  # Node does not yet exist
            if i == length(kmer)  # Node is a leaf, need to add a count entry
                push!(trie.counts, 0)
                # println(trie.counts)
                next = Node{Sb}(n, Node{Sb}[], length(trie.counts))
            else
                next = Node{Sb}(n, Node{Sb}[], 0)
            end
            push!(nexts, next)
        end
        last = next
        nexts = next.nexts
    end
    trie.counts[last.count_idx] += 1
    
end

function get_next(v::Vector{Node{Sb}}, symbol::Sb) where {Sb<:BioSymbol}
    for n in v
        n.value == symbol && return n
    end
end

function init_trie(K::Int, Sb::Type{<:BioSymbol}=DNA)
    alpha = Node{Sb}[Node{Sb}(symbol, Node{Sb}[], 0) for symbol in alphabet(Sb) if !(isambiguous(symbol) || isgap(symbol))]
    return Trie{Sb}(K, alpha, UInt32[])
end

function twobit(n::NucleicAcid)::UInt8
    bits = compatbits(n)
    return UInt8(trailing_zeros(bits))
end

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

# function Kmers.build_kmer(
#     ::Union{Base.HasLength, Base.HasShape},
#     ::Kmers.RecodingScheme,
#     T::Type,
#     s::Vector{AminoAcid},
# )
#     # println("WE GOT IT")
#     length(s) == Kmers.ksize(T) || error("Length of sequence must be K elements to build Kmer")
#     data = Kmers.zero_tuple(T)
#     A = Alphabet(T)
#     bps = BioSequences.bits_per_symbol(A)
#     for element in s
#         symbol = convert(eltype(A), element)
#         isnothing(BioSequences.encode(A, symbol)) && println(s)
#         carry = UInt(BioSequences.encode(A, symbol))
#         (_, data) = Kmers.leftshift_carry(data, bps, carry)
#     end
#     T(Kmers.unsafe, data)
# end

function k_merize(sequence::LongSequence{Ab}; K::Int) where {Ab<:Alphabet}
    length(sequence) < K && error("Read of size $(length(sequence)) too small for $K-mers")
    to_return = Kmer{Ab, K}[]
    for i in 1:(length(sequence)-K+1)
        push!(to_return, Kmer{Ab, K}(sequence[i:i+K-1]))
    end
    return to_return
end

function k_merize(sequence::LongAA; K::Int)
    length(sequence) < K && error("Read of size $(length(sequence)) too small for $K-mers")
    to_return = Kmer{AAAlphabet, K}[]
    for i in 1:(length(sequence)-K+1)
        push!(to_return, Kmer{AAAlphabet, K}(sequence[i:i+K-1]))
    end
    return to_return
end

function jello_hash(fastq::String; K::Int)
    hash = Dict{UInt64, UInt32}()
    wcl = countlines(fastq)
    progress = Progress(wcl, desc="Parsing fastq into hash...")
    open(fastq) do f
        for (i, l) in enumerate(eachline(f))
            (i-2)%4 != 0 && continue    # Skips non-sequence lines of fastq
            
            'N' in l && continue
            seq = LongSequence{DNAAlphabet{2}}(l)
            kmers = translate.(k_merize(seq, K=K))
            for kmer in kmers
                @time AminoAcid('*') in kmer && continue
                if !(haskey(hash, kmer.data[1]))
                    hash[kmer.data[1]] = UInt32(0)
                end
                hash[kmer.data[1]] += 1
            end
            next!(progress; showvalues=[
                ("unique k-mers", length(hash))])
        end
    end
    return hash
end

function Kct(hash_table::Dict{UInt64, UInt32}, K::Int, A::Ab=DNAAlphabet{2}) where {Ab<:Type{<:Alphabet}}

    count_type = UInt32
    symbol_size = BioSequences.bits_per_symbol(A())
    required_chunks = ceil(Int64, K/(64/symbol_size))  # On 64 bits

    nb = Vector{count_type}()
    kct = Kct(1, K)

    @showprogress desc="parsing hash table into KCT..." for (k_bits, k_val) in hash_table
        tmp_seq = Kmer{A, K, required_chunks}(Kmers.unsafe, (k_bits,))
        push!(kct.table, KRecord(tmp_seq, UInt32(0)))
        push!(nb, k_val)
    end
    
    _dedup(kct, nb)
    sort!(kct)
    return kct
end

# abstract type Extension{E} end
abstract type Fq end

struct FastqGZ <: Fq
    iterable::GZipStream
end

struct Fastq <: Fq
    iterable::IOStream
end

struct Mzid
    iterable::Vector{EzXML.Node}
end

function Base.iterate(seqfile::F, state=nothing) where {F<:Fq}
    eof(seqfile.iterable) && return nothing
    readline(seqfile.iterable)
    seq = strip(readline(seqfile.iterable))
    readline(seqfile.iterable)
    readline(seqfile.iterable)
    return seq, state
end

function Base.iterate(seqfile::Mzid, state::Int=1)
    return state < length(seqfile.iterable) ? (seqfile.iterable[state].content, state+1) : nothing
end

stream(::Type{Val{:fastq}}, path::String) = Fastq(open(path, "r"))
stream(::Type{Val{:gz}}, path::String) = FastqGZ(gzopen(path, "r"))  # TODO: Assumes fastq.gz
stream(::Type{Val{:mzid}}, path::String) = Mzid(findall("//*[local-name()='Seq']", readxml(path)))
stream(path::String) = stream(Val{Symbol(split(path, ".")[end])}, path)


"""
    chunk_stream(file::String, K::Int, merge_queue::Channel{Dict{UInt64, UInt32}}, chunking::Int=1_000_000)

This function does a thing
"""
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
            push!(counting_tasks, @spawn process_thread(to_process, K, merge_queue,  verbose=true))
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
        push!(counting_tasks, @spawn process_thread(current_chunk, K, merge_queue, verbose=true))
        !verbose && ProgressMeter.update!(progress; showvalues=[
            ("active counting tasks", length(filter(x->x.state!=:done, counting_tasks))),
            ("items in merge queue", merge_queue.n_avail_items),
            ("chunks sent", i÷chunking+1)
        ])
    end
    wait.(counting_tasks)
    finish!(progress)
end

### DEPRECATED
function jello_threaded_hash(fastq::String, K::Int, chunking::Int=1_000_000, queue_size::Int=32)

    merge_queue = Channel{Dict{UInt64, UInt32}}(queue_size)
    mother_hash = Dict{UInt64, UInt32}()
    merger = @async try
        for child_hash in merge_queue
            mother_hash = merge(+, mother_hash, child_hash)
        end
    catch e
        @error "Merger task crashed" exception=(e, catch_backtrace())
    end

    chunk_stream(fastq, K, merge_queue, chunking)

    close(merge_queue)
    wait(merger)
    return mother_hash
end

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

function process_thread(chunk::Vector{String}, K::Int, merge_queue::Channel{Dict{UInt64, UInt32}}; verbose::Bool=false)
    hash = Dict{UInt64, UInt32}()
    start = now()
    try 
        for l in chunk
            'N' in l && continue
            seq = LongSequence{DNAAlphabet{2}}(l)
            length(seq) < K && continue
            kmers = translate.(k_merize(seq, K=K))
            for kmer in kmers
                isnothing(kmer) && continue
                # AminoAcid('*') in kmer && continue
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
    # verbose && tprintln("Completed counting in $(now()-start)")
end

### EXPERIMENTAL

function jello_trie(fastq::String, Sb::Type{<:BioSymbol}=DNA; K::Int)
    # GC.enable(false)
    trie = init_trie(K, Sb)
    wcl = countlines(fastq)
    progress = Progress(wcl, desc="Parsing fastq into trie...")
    open(fastq) do f
        for (i, l) in enumerate(eachline(f))
            (i-2)%4 != 0 && continue    # Skips non-sequence lines of fastq
            
            seq = bioseq(l)
            kmers = k_merize(seq, K=K)
            push!.(Ref(trie), kmers)
            next!(progress; showvalues=[
                ("unique k-mers", length(trie.counts))
                ("counts memory size", Base.format_bytes(Base.summarysize(trie.counts)))
                ])
        end
    end
    # GC.enable(true)
    return trie
end

### DEPRECATED

# function jello_kct_chunk(fastq::String; K::Int, chunking::Int=1_000_000, queue_size::Int=32)
#     lines = readlines(open(fastq, "r"))[2:4:end]
#     chunks = [lines[i:min(i+chunking-1, end)] for i in 1:chunking:length(lines)]
#     merge_queue = Channel{Kct{1, K}}(queue_size)
#     mother_kct = Kct(1, K)
#     # Merger thread. Works asynchronously from workers.
#     merger = @async try
#         for child_kct in merge_queue
#             mother_kct = merge(mother_kct, child_kct)
#         end
#     catch e
#         @error "Merger task crashed" exception=(e, catch_backtrace())
#     end

#     progress = Progress(length(chunks), desc = "Processing $chunking k-mers chunks...")
#     # Worker threads, building child Kcts asynchronously
#     @threads for chunk in chunks
#         child_kct = Kct(1, K)
#         for l in chunk
#             # println(l)
#             seq = bioseq(l)
#             DNA_N in seq && continue
#             k_mers = k_merize(seq, K=K)
#             for k_mer in k_mers
#                 buf = JaggedRLEArray(UInt32)
#                 push!(buf, UInt32[1])
#                 id = push!(child_kct, buf)
#                 # println(k_mer)
#                 push!(child_kct.table, KRecord(k_mer, id))
#             end
#         end
#         # println("starting a sort")
#         sort!(child_kct)
#         # println(merge_queue.n_avail_items)
#         put!(merge_queue, child_kct)
#         # println("finished one chunk")
#         next!(progress; showvalues=[
#                 ("items in merge queue", merge_queue.n_avail_items)
#                 ])
#     end

#     close(merge_queue)
#     wait(merger)

#     return mother_kct
# end