# Global loader. Should work regardless of kct version
# Dispatches via version, sorted in the first 64 bits of kct file.
function load_kct(path::String)
    open(path, "r") do io
        version = read(io, Float64)
        kct = load(io, Val(version))
        return kct
    end
end

# Global writer for KCTs. Dispatches the writing strategy based on kct.version
# Version is always written in the first 64 bits (Float64) to dispatch loader
function write_kct(kct::NeoKCT{K, Ab}, path::String) where {K, Ab}
    open(path, "w") do io
        write(io, kct.version)
        write(io, kct, Val(kct.version))
    end
end

### VERSION DEPENDENT LOADERS ###

## V1.1 ##
function Base.write(io::IO, kct::NeoKCT{K, Ab}, version::V) where {K, Ab<:Alphabet, V<:Val{1.1}}
    Ab_name = String(Ab.name.singletonname)
    kmer_C = isempty(kct.table) ? 1 : length(kct.table[1].seq.data)
    write(io, Int64(K))
    write(io, Int64(length(Ab_name)))
    write(io, codeunits(Ab_name))
    write(io, Int64(kmer_C))
    write(io, Int64(length(kct.table)))
    write(io, Int64(length(kct.counts.words)))
    write(io, Int64(length(kct.counts.bitmap)))
    write(io, Int64(kct.samples.x))
    for k_elem in kct.table
        for d in k_elem.seq.data; write(io, d); end
        write(io, Int32(length(k_elem.chunk_ids)))
        for cid in k_elem.chunk_ids; write(io, cid); end
    end
    write(io, kct.counts.words)  # bulk write
    for chunk in kct.counts.bitmap.chunks
        write(io, chunk)
    end
end

function load(io::IO, version::Val{1.1})
    K = read(io, Int64)
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))  # once only
    kmer_C = read(io, Int64)
    table_length = read(io, Int64)
    words_length = read(io, Int64)
    bitmap_length = read(io, Int64)
    samples_count = read(io, Int64)

    table = Vector{K_Element{K, Ab, kmer_C}}(undef, table_length)
    @showprogress "Reading K-mer Table..." for i in 1:table_length
        data  = NTuple{kmer_C, UInt64}(read(io, UInt64) for _ in 1:kmer_C)
        word_C = Int64(read(io, Int32))
        chunk_ids = NTuple{word_C, UInt32}(read(io, UInt32) for _ in 1:word_C)
        table[i] = K_Element{K, Ab, kmer_C}(Kmer{Ab, K, kmer_C}(Kmers.unsafe, data), chunk_ids)
    end

    words = Vector{UInt64}(undef, words_length)
    read!(io, words)  # bulk read

    bitmap = BitVector(undef, bitmap_length)
    for i in eachindex(bitmap.chunks)
        bitmap.chunks[i] = read(io, UInt64)
    end

    idx = Ref(0)=>fill(0:-1, 4^15)
    kct = NeoKCT(table, PackedArray{UInt32, UInt64}(words, bitmap, words_length), idx, Ref(samples_count), 1.1)
    compute_index!(kct)
    return kct
end
## V1.0 ##
function Base.write(io::IO, tuple::Tuple, version::Union{Val{1.0}, Val{1.1}})
    write(io, length(tuple))
    write.(io, tuple)
end

function Base.write(io::IO, k_elem::K_Element{K, Ab, C}, version::Union{Val{1.0}, Val{1.1}}) where {K, Ab<:AAAlphabet, C}
    write(io, K)
    write(io, codeunits(String(Ab.name.singletonname)))
    write(io, C)
    write.(io, k_elem.seq.data)
    write(io, k_elem.chunk_ids, version)
end

function Base.write(io::IO, kct::NeoKCT{K, Ab}, version::V) where {K, Ab<:Alphabet, V<:Val{1.0}}
    write(io, kct.version)
    write(io, K)
    write(io, (String(Ab.name.singletonname)))
    write(io, length(kct.table))
    write(io, length(kct.counts.words))
    write(io, length(kct.counts.bitmap))
    write.(io, kct.table, version)
    write(io, kct.counts.words)
    write(io, kct.counts.bitmap)
    write(io, kct.samples.x)
end

function load(io::IO, version::Val{1.0})
    K = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:10])))  # TODO: Fix this abomination
    table_length = read(io, Int64)
    words_length = read(io, Int64)
    bitmap_length = read(io, Int64)
    table = Array{K_Element{K, Ab, 1}, 1}(undef, table_length)
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
    bitmap = BitVector(undef, bitmap_length)
    for i in eachindex(bitmap.chunks)
        bitmap.chunks[i] = read(io, UInt64)
    end
    samples = Ref(read(io, Int64))
    idx = Ref(20)=>fill(0:-1, 4^15)
    return NeoKCT(table, PackedArray{UInt32, UInt64}(words, bitmap, words_length), idx, samples, 1.0)
end