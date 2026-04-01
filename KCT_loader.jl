function load_kct(path::String)
    open(path, "r") do io
        version = read(io, Float64)
        kct = load(io, version=version)
        return kct
    end
end

function write_kct(kct::NeoKCT{K, Ab}, path::String)
    open(path, "w") do io
        write(io, kct, version = kct.version)
    end
end

### VERSION DEPENDENT LOADERS ###

## V1.0 ##
function Base.write(io::IO, tuple::Tuple; version::V) where {V in [1.0]}
    write(io, length(tuple))
    write.(io, tuple)
end

function Base.write(io::IO, k_elem::K_Element{K, Ab, C}; version::V) where {K, Ab<:AAAlphabet, C, V in [1.0]}
    write(io, K)
    write(io, codeunits(String(Ab.name.singletonname)))
    write(io, C)
    write.(io, k_elem.seq.data, version=version)
    write(io, k_elem.chunk_ids)
end

function Base.write(io::IO, kct::NeoKCT{K, Ab}; version::V) where {K, Ab<:Alphabet, V in [1.0]}
    write(io, kct.version)
    write(io, K)
    write(io, (String(Ab.name.singletonname)))
    write(io, length(kct.table))
    write(io, length(kct.counts.words))
    write(io, length(kct.counts.bitmap))
    write.(io, kct.table, version=version)
    write(io, kct.counts.words)
    write(io, kct.counts.bitmap)
    write(io, kct.samples.x)
end

function load(io::IO; version::V) where {V in [1.0]}
    io = open(path, "r")
    version = read(io, Float64)
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
    close(io)
    idx = fill(0:-1, 4^15)
    return NeoKCT(table, PackedArray{UInt32, UInt64}(words, bitmap, words_length), idx, samples, version)

end