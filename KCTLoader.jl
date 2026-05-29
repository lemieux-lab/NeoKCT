# Global loader. Dispatches via version stored in the first 64 bits of the file.
function load_kct(path::String)
    open(path, "r") do io
        version = read(io, Float64)
        return load(io, Val(version))
    end
end

# Global writer. Dispatches the writing strategy based on kct.version.
function write_kct(kct::NeoKCT{K, Ab}, path::String) where {K, Ab}
    open(path, "w") do io
        write(io, kct.version)
        write(io, kct, Val(kct.version))
    end
end

function write_kct(rich::RichKCT{K, Ab}, path::String) where {K, Ab}
    open(path, "w") do io
        write(io, rich.version)
        write(io, rich, Val(rich.version))
    end
end

const _WORD_TYPES = Dict{Int64, DataType}(
    1 => UInt8, 2 => UInt16, 4 => UInt32, 8 => UInt64, 16 => UInt128
)

function get_version(path::String)
    open(path, "r") do io
        return read(io, Float64)
    end
end

### GENOMIC INDEX SERIALISATION ###

const GIDX_VERSION = 1.0

function write_gidx(gidx::GenomicIndex{K, Ab}, path::String) where {K, Ab}
    open(path, "w") do io
        write(io, GIDX_VERSION)
        _write_gidx(io, gidx, Val(GIDX_VERSION))
    end
end

function load_gidx(path::String)
    open(path, "r") do io
        version = read(io, Float64)
        return _load_gidx(io, Val(version))
    end
end

function _write_gidx(io::IO, gidx::GenomicIndex{K, Ab}, ::Val{1.0}) where {K, Ab}
    Ab_name = String(Ab.name.singletonname)
    n_kmers = length(gidx.seqs)
    write(io, Int64(K))
    write(io, Int64(length(Ab_name))); write(io, codeunits(Ab_name))
    write(io, Int64(n_kmers))
    write(io, Int64(sizeof(eltype(gidx.seqs.checkpoints))))
    write(io, Int64(sizeof(eltype(gidx.seqs.deltas))))
    write(io, Int64(gidx.seqs.checkpoint_interval))
    write(io, Int64(length(gidx.seqs.checkpoints)))
    write(io, Int64(length(gidx.seqs.regular_cp_idx)))
    write(io, gidx.seqs.checkpoints)
    write(io, gidx.seqs.deltas)
    write(io, Int64.(gidx.seqs.regular_cp_idx))
    write(io, Int64(length(gidx.biotype_names)))
    for name in gidx.biotype_names
        write(io, Int64(length(name))); write(io, codeunits(name))
    end
    write(io, Int64(length(gidx.pool)))
    write(io, gidx.pool)
    write(io, gidx.biotype_ids)
end

function _load_gidx(io::IO, ::Val{1.0})
    K = Int(read(io, Int64))
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    n_kmers = read(io, Int64)
    C_type = _WORD_TYPES[read(io, Int64)]
    D_type = _WORD_TYPES[read(io, Int64)]
    cp_interval = Int(read(io, Int64))
    n_cp = read(io, Int64)
    n_rci = read(io, Int64)
    checkpoints = Vector{C_type}(undef, n_cp); read!(io, checkpoints)
    deltas = Vector{D_type}(undef, n_kmers); read!(io, deltas)
    rci_i64 = Vector{Int64}(undef, n_rci); read!(io, rci_i64)
    seqs = DeltaArray(checkpoints, deltas, UInt64.(rci_i64), cp_interval)
    n_names = read(io, Int64)
    biotype_names = Vector{String}(undef, n_names)
    for i in 1:n_names
        len = read(io, Int64)
        biotype_names[i] = String([read(io, UInt8) for _ in 1:len])
    end
    pool_len = read(io, Int64)
    pool = Vector{UInt64}(undef, pool_len); read!(io, pool)
    biotype_ids = Vector{UInt16}(undef, n_kmers); read!(io, biotype_ids)
    gidx = GenomicIndex{K, Ab}(seqs, biotype_ids, pool, Ref(20) => fill(0:-1, 4^15), biotype_names)
    compute_index!(gidx)
    return gidx
end

### VERSION DEPENDENT LOADERS ###
## V2.0 — RichKCT (NeoKCT v1.4 body + BiotypLayer) ##
function Base.write(io::IO, rich::RichKCT{K, Ab, C}, ::Val{2.0}) where {K, Ab<:Alphabet, C}
    write(io, rich.kct, Val(1.4))  # NeoKCT body, no version prefix
    write(io, Int64(length(rich.biotypes.biotype_names)))
    for name in rich.biotypes.biotype_names
        write(io, Int64(length(name))); write(io, codeunits(name))
    end
    write(io, Int64(length(rich.biotypes.pool)))
    write(io, rich.biotypes.pool)
    write(io, rich.biotypes.ids)
end

function load(io::IO, ::Val{2.0})
    kct = load(io, Val(1.4))
    n_names = read(io, Int64)
    names = Vector{String}(undef, n_names)
    for i in 1:n_names
        len = read(io, Int64)
        names[i] = String([read(io, UInt8) for _ in 1:len])
    end
    pool_len = read(io, Int64)
    pool = Vector{UInt64}(undef, pool_len); read!(io, pool)
    ids = Vector{UInt16}(undef, length(kct.seqs)); read!(io, ids)
    return _wrap_rich(kct, BiotypLayer(ids, pool, names))
end

# Dispatch helper to extract {K,Ab,C} from the loaded NeoKCT.
_wrap_rich(kct::NeoKCT{K,Ab,C}, biotypes::BiotypLayer) where {K,Ab,C} =
    RichKCT{K,Ab,C}(kct, biotypes)

## V1.4 ##
function Base.write(io::IO, kct::NeoKCT{K, Ab, C}, ::Val{1.4}) where {K, Ab<:Alphabet, C}
    W = eltype(kct.counts.words)
    Ab_name = String(Ab.name.singletonname)
    write(io, Int64(K))
    write(io, Int64(length(Ab_name))); write(io, codeunits(Ab_name))
    write(io, Int64(length(kct.seqs)))  # n_kmers
    write(io, Int64(length(kct.flat_cids)))
    write(io, Int64(length(kct.counts.words)))
    write(io, Int64(length(kct.counts.bitmap)))
    write(io, Int64(kct.samples.x))
    write(io, Int64(sizeof(W)))
    # DeltaArray fields
    write(io, Int64(kct.seqs.checkpoint_interval))
    write(io, Int64(length(kct.seqs.checkpoints)))
    write(io, Int64(length(kct.seqs.regular_cp_idx)))
    write(io, kct.seqs.checkpoints)
    write(io, kct.seqs.deltas)  # Vector{UInt32}
    write(io, Int64.(kct.seqs.regular_cp_idx))
    # rest
    write(io, kct.n_cids)
    write(io, kct.flat_cids)
    write(io, kct.counts.words)
    for chunk in kct.counts.bitmap.chunks; write(io, chunk); end
end

function load(io::IO, ::Val{1.4})
    K = read(io, Int64)
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    n_kmers = read(io, Int64)
    n_flat_cids = read(io, Int64)
    words_len = read(io, Int64)
    bitmap_len = read(io, Int64)
    n_samples = read(io, Int64)
    word_size_b = read(io, Int64)
    W = _WORD_TYPES[word_size_b]
    # DeltaArray fields
    C = Int(read(io, Int64))
    n_cp = read(io, Int64)
    n_rci = read(io, Int64)
    checkpoints = Vector{UInt64}(undef, n_cp); read!(io, checkpoints)
    deltas = Vector{UInt32}(undef, n_kmers); read!(io, deltas)
    rci_i64 = Vector{Int64}(undef, n_rci); read!(io, rci_i64)
    regular_cp_idx = Vector{Int}(rci_i64)
    regular_cp_idx = UInt64.(regular_cp_idx)  # TODO: Fix Int → Uint Read/Write (now UInt for regular_cp_index)
    # rest
    n_cids = Vector{UInt16}(undef, n_kmers); read!(io, n_cids)
    flat_cids = Vector{UInt32}(undef, n_flat_cids); read!(io, flat_cids)
    words = Vector{W}(undef, words_len); read!(io, words)
    bitmap = BitVector(undef, bitmap_len)
    for i in eachindex(bitmap.chunks); bitmap.chunks[i] = read(io, UInt64); end

    seqs = DeltaArray(checkpoints, deltas, regular_cp_idx, C)
    idx = Ref(0) => fill(0:-1, 4^15)
    kct = NeoKCT{K, Ab, 1}(seqs, n_cids, flat_cids,
               PackedArray{UInt32, W}(words, bitmap),
               idx, Ref(n_samples), 1.4)
    compute_index!(kct)
    return kct
end

## V1.3 ##
function Base.write(io::IO, kct::NeoKCT{K, Ab, C}, ::Val{1.3}) where {K, Ab<:Alphabet, C}
    W = eltype(kct.counts.words)
    seqs = collect(kct.seqs)   # decode to flat for legacy format
    Ab_name = String(Ab.name.singletonname)
    write(io, Int64(K))
    write(io, Int64(length(Ab_name)));  write(io, codeunits(Ab_name))
    write(io, Int64(length(seqs)))
    write(io, Int64(length(kct.flat_cids)))
    write(io, Int64(length(kct.counts.words)))
    write(io, Int64(length(kct.counts.bitmap)))
    write(io, Int64(kct.samples.x))
    write(io, Int64(sizeof(W)))
    write(io, seqs)
    write(io, kct.n_cids)
    write(io, kct.flat_cids)
    write(io, kct.counts.words)
    for chunk in kct.counts.bitmap.chunks; write(io, chunk); end
end

function load(io::IO, ::Val{1.3})
    K = read(io, Int64)
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    n_kmers = read(io, Int64)
    n_flat_cids = read(io, Int64)
    words_len = read(io, Int64)
    bitmap_len = read(io, Int64)
    n_samples = read(io, Int64)
    word_size_b = read(io, Int64)
    W = _WORD_TYPES[word_size_b]

    raw_seqs = Vector{UInt64}(undef, n_kmers); read!(io, raw_seqs)
    n_cids = Vector{UInt16}(undef, n_kmers); read!(io, n_cids)
    flat_cids = Vector{UInt32}(undef, n_flat_cids); read!(io, flat_cids)
    words = Vector{W}(undef, words_len); read!(io, words)
    bitmap = BitVector(undef, bitmap_len)
    for i in eachindex(bitmap.chunks); bitmap.chunks[i] = read(io, UInt64); end

    seqs = DeltaArray(raw_seqs, DEFAULT_CHECKPOINT_INTERVAL)
    idx = Ref(0) => fill(0:-1, 4^15)
    kct = NeoKCT{K, Ab, 1}(seqs, n_cids, flat_cids,
               PackedArray{UInt32, W}(words, bitmap),
               idx, Ref(n_samples), 1.4)
    compute_index!(kct)
    return kct
end

## V1.2 ##
function load(io::IO, ::Val{1.2})
    K = read(io, Int64)
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    n_kmers = read(io, Int64)
    n_flat_cids = read(io, Int64)
    words_len = read(io, Int64)
    bitmap_len = read(io, Int64)
    n_samples = read(io, Int64)
    word_size_b = read(io, Int64)
    W = _WORD_TYPES[word_size_b]

    raw_seqs = Vector{UInt64}(undef, n_kmers); read!(io, raw_seqs)
    offsets_u32 = Vector{UInt32}(undef, n_kmers + 1); read!(io, offsets_u32)
    n_cids = UInt16.(diff(offsets_u32))
    flat_cids = Vector{UInt32}(undef, n_flat_cids); read!(io, flat_cids)
    words = Vector{W}(undef, words_len); read!(io, words)
    bitmap = BitVector(undef, bitmap_len)
    for i in eachindex(bitmap.chunks); bitmap.chunks[i] = read(io, UInt64); end

    seqs = DeltaArray(raw_seqs, DEFAULT_CHECKPOINT_INTERVAL)
    idx = Ref(0) => fill(0:-1, 4^15)
    kct = NeoKCT{K, Ab, 1}(seqs, n_cids, flat_cids,
               PackedArray{UInt32, W}(words, bitmap),
               idx, Ref(n_samples), 1.4)
    compute_index!(kct)
    return kct
end

# ## V1.1 ##  /!\ NOW FULLY DEPRECATED /!\
# function Base.write(io::IO, kct::NeoKCT{K, Ab}, version::V) where {K, Ab<:Alphabet, V<:Val{1.1}}
#     Ab_name = String(Ab.name.singletonname)
#     kmer_C = isempty(kct.table) ? 1 : length(kct.table[1].seq.data)
#     write(io, Int64(K))
#     write(io, Int64(length(Ab_name)))
#     write(io, codeunits(Ab_name))
#     write(io, Int64(kmer_C))
#     write(io, Int64(length(kct.table)))
#     write(io, Int64(length(kct.counts.words)))
#     write(io, Int64(length(kct.counts.bitmap)))
#     write(io, Int64(kct.samples.x))
#     for k_elem in kct.table
#         for d in k_elem.seq.data; write(io, d); end
#         write(io, Int32(length(k_elem.chunk_ids)))
#         for cid in k_elem.chunk_ids; write(io, cid); end
#     end
#     write(io, kct.counts.words)  # bulk write
#     for chunk in kct.counts.bitmap.chunks
#         write(io, chunk)
#     end
# end

# function load(io::IO, version::Val{1.1})
#     K = read(io, Int64)
#     Ab_name_len = read(io, Int64)
#     Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))  # once only
#     kmer_C = read(io, Int64)
#     table_length = read(io, Int64)
#     words_length = read(io, Int64)
#     bitmap_length = read(io, Int64)
#     samples_count = read(io, Int64)

#     table = Vector{K_Element{K, Ab, kmer_C}}(undef, table_length)
#     @showprogress "Reading K-mer Table..." for i in 1:table_length
#         data  = NTuple{kmer_C, UInt64}(read(io, UInt64) for _ in 1:kmer_C)
#         word_C = Int64(read(io, Int32))
#         chunk_ids = NTuple{word_C, UInt32}(read(io, UInt32) for _ in 1:word_C)
#         table[i] = K_Element{K, Ab, kmer_C}(Kmer{Ab, K, kmer_C}(Kmers.unsafe, data), chunk_ids)
#     end

#     words = Vector{UInt64}(undef, words_length)
#     read!(io, words)  # bulk read

#     bitmap = BitVector(undef, bitmap_length)
#     for i in eachindex(bitmap.chunks)
#         bitmap.chunks[i] = read(io, UInt64)
#     end

#     idx = Ref(0)=>fill(0:-1, 4^15)
#     kct = NeoKCT(table, PackedArray{UInt32, UInt64}(words, bitmap, words_length), idx, Ref(samples_count), 1.1)
#     compute_index!(kct)
#     return kct
# end
# ## V1.0 ##
# function Base.write(io::IO, tuple::Tuple, version::Union{Val{1.0}, Val{1.1}})
#     write(io, length(tuple))
#     write.(io, tuple)
# end

# function Base.write(io::IO, k_elem::K_Element{K, Ab, C}, version::Union{Val{1.0}, Val{1.1}}) where {K, Ab<:AAAlphabet, C}
#     write(io, K)
#     write(io, codeunits(String(Ab.name.singletonname)))
#     write(io, C)
#     write.(io, k_elem.seq.data)
#     write(io, k_elem.chunk_ids, version)
# end

# function Base.write(io::IO, kct::NeoKCT{K, Ab}, version::V) where {K, Ab<:Alphabet, V<:Val{1.0}}
#     write(io, kct.version)
#     write(io, K)
#     write(io, (String(Ab.name.singletonname)))
#     write(io, length(kct.table))
#     write(io, length(kct.counts.words))
#     write(io, length(kct.counts.bitmap))
#     write.(io, kct.table, version)
#     write(io, kct.counts.words)
#     write(io, kct.counts.bitmap)
#     write(io, kct.samples.x)
# end

# function load(io::IO, version::Val{1.0})
#     K = read(io, Int64)
#     Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:10])))  # TODO: Fix this abomination
#     table_length = read(io, Int64)
#     words_length = read(io, Int64)
#     bitmap_length = read(io, Int64)
#     table = Array{K_Element{K, Ab, 1}, 1}(undef, table_length)
#     @showprogress for i in 1:table_length
#         K = read(io, Int64)
#         Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:10])))  # TODO: Fix this abomination
#         kmer_C = read(io, Int64)
#         data = Tuple(read(io, UInt64) for _ in 1:kmer_C)
#         word_C = read(io, Int64)
#         chunk_ids = Tuple(read(io, UInt32) for _ in 1:word_C)
#         table[i] = K_Element{K, Ab, kmer_C}(Kmer{Ab, K, kmer_C}(Kmers.unsafe, data), chunk_ids)
#     end
#     words = Vector{UInt64}(undef, words_length)
#     @showprogress for i in 1:words_length
#         words[i] = read(io, UInt64)
#     end
#     bitmap = BitVector(undef, bitmap_length)
#     for i in eachindex(bitmap.chunks)
#         bitmap.chunks[i] = read(io, UInt64)
#     end
#     samples = Ref(read(io, Int64))
#     idx = Ref(20)=>fill(0:-1, 4^15)
#     return NeoKCT(table, PackedArray{UInt32, UInt64}(words, bitmap, words_length), idx, samples, 1.0)
# end