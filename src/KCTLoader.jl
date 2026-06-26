const KCT_VERSION = 3.0

const _WORD_TYPES = Dict{Int64, DataType}(
    1 => UInt8, 2 => UInt16, 4 => UInt32, 8 => UInt64, 16 => UInt128
)

## PUBLIC API ##

function write_kct(kct::KCT{K, Ab}, path::String) where {K, Ab}
    open(path, "w") do io
        write(io, KCT_VERSION)
        _write_kct(io, kct, Val(KCT_VERSION))
    end
end

function load_kct(path::String)
    open(path, "r") do io
        version = read(io, Float64)
        return _load_kct(io, Val(version))
    end
end

function get_version(path::String)
    open(path, "r") do io
        return read(io, Float64)
    end
end

## V3.0 FORMAT ##
# Header: [Float64 version][Int64 K][Int64 Ab_name_len][UInt8... Ab_name][UInt8 layers_mask][Int64 n_kmers]
# KmerLayer: [Int64 sizeof(C)][Int64 sizeof(D)][Int64 cp_interval][Int64 n_cp][Int64 n_rci][C... cps][D... deltas][Int64... rci]
# CountsLayer (mask bit 0): [Int64 n_samples][Int64 n_flat_cids][Int64 words_len][Int64 bitmap_len][Int64 sizeof(W)][UInt16... n_cids][UInt32... flat_cids][W... words][UInt64... bitmap.chunks]
# BiotypLayer (mask bit 1): [Int64 n_names]([Int64 len][UInt8... name]...)[Int64 pool_len][UInt64... pool][UInt16... ids]

_layers_mask(::Nothing, ::Nothing) = UInt8(0)
_layers_mask(::CountsLayer, ::Nothing) = UInt8(1)
_layers_mask(::Nothing, ::BiotypLayer) = UInt8(2)
_layers_mask(::CountsLayer, ::BiotypLayer) = UInt8(3)

function _write_kct(io::IO, kct::KCT{K, Ab, Counts, Biotype, C, D}, ::Val{3.0}) where {K, Ab<:Alphabet, Counts, Biotype, C<:Unsigned, D<:Unsigned}
    Ab_name = String(Ab.name.singletonname)
    write(io, Int64(K))
    write(io, Int64(length(Ab_name))); write(io, codeunits(Ab_name))
    write(io, _layers_mask(kct.counts, kct.biotype))
    write(io, Int64(length(kct.kmer.seqs)))
    write(io, Int64(sizeof(C)))
    write(io, Int64(sizeof(D)))
    write(io, Int64(kct.kmer.seqs.checkpoint_interval))
    write(io, Int64(length(kct.kmer.seqs.checkpoints)))
    write(io, Int64(length(kct.kmer.seqs.regular_cp_idx)))
    write(io, kct.kmer.seqs.checkpoints)
    write(io, kct.kmer.seqs.deltas)
    write(io, Int64.(kct.kmer.seqs.regular_cp_idx))
    _write_counts(io, kct.counts)
    _write_biotype(io, kct.biotype)
end

_write_counts(::IO, ::Nothing) = nothing
function _write_counts(io::IO, cl::CountsLayer)
    W = eltype(cl.counts.words)
    write(io, Int64(cl.samples.x))
    write(io, Int64(length(cl.flat_cids)))
    write(io, Int64(length(cl.counts.words)))
    write(io, Int64(length(cl.counts.bitmap)))
    write(io, Int64(sizeof(W)))
    write(io, cl.n_cids)
    write(io, cl.flat_cids)
    write(io, cl.counts.words)
    for chunk in cl.counts.bitmap.chunks; write(io, chunk); end
end

_write_biotype(::IO, ::Nothing) = nothing
function _write_biotype(io::IO, bl::BiotypLayer)
    write(io, Int64(length(bl.biotype_names)))
    for name in bl.biotype_names
        write(io, Int64(length(name))); write(io, codeunits(name))
    end
    write(io, Int64(length(bl.pool)))
    write(io, bl.pool)
    write(io, bl.ids)
end

function _load_kct(io::IO, ::Val{3.0})
    K = Int(read(io, Int64))
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    layers_mask = read(io, UInt8)
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
    kl = KmerLayer{K, Ab, C_type, D_type}(seqs, Ref(20) => fill(0:-1, 1 << max(0, (K - 4) * bits_per_symbol(Ab()))))
    cl = (layers_mask & UInt8(1)) != 0 ? _read_counts(io, n_kmers) : nothing
    bl = (layers_mask & UInt8(2)) != 0 ? _read_biotype(io, n_kmers) : nothing
    kct = _build_kct(kl, cl, bl)
    compute_index!(kct)
    return kct
end

function _read_counts(io::IO, n_kmers::Int64)
    n_samples = read(io, Int64)
    n_flat_cids = read(io, Int64)
    words_len = read(io, Int64)
    bitmap_len = read(io, Int64)
    W = _WORD_TYPES[read(io, Int64)]
    n_cids = Vector{UInt16}(undef, n_kmers); read!(io, n_cids)
    flat_cids = Vector{UInt32}(undef, n_flat_cids); read!(io, flat_cids)
    words = Vector{W}(undef, words_len); read!(io, words)
    bitmap = BitVector(undef, bitmap_len)
    for i in eachindex(bitmap.chunks); bitmap.chunks[i] = read(io, UInt64); end
    return CountsLayer(n_cids, flat_cids, PackedArray{UInt32, W}(words, bitmap), Ref(n_samples))
end

function _read_biotype(io::IO, n_kmers::Int64)
    n_names = read(io, Int64)
    biotype_names = Vector{String}(undef, n_names)
    for i in 1:n_names
        len = read(io, Int64)
        biotype_names[i] = String([read(io, UInt8) for _ in 1:len])
    end
    pool_len = read(io, Int64)
    pool = Vector{UInt64}(undef, pool_len); read!(io, pool)
    ids = Vector{UInt16}(undef, n_kmers); read!(io, ids)
    return BiotypLayer(ids, pool, biotype_names)
end

_build_kct(kl, ::Nothing, ::Nothing) = KCT(kl)
_build_kct(kl, cl::CountsLayer, ::Nothing) = KCT(kl, cl)
_build_kct(kl, cl::CountsLayer, bl::BiotypLayer) = KCT(kl, cl, bl)
_build_kct(kl, ::Nothing, bl::BiotypLayer) = KCT(kl, bl)

## RETROCOMPAT LOADERS ##

# V2.0 (RichKCT) → KCT{K,Ab,CountsLayer,BiotypLayer}
function _load_kct(io::IO, ::Val{2.0})
    kct = _load_kct(io, Val(1.4))
    n_names = read(io, Int64)
    names = Vector{String}(undef, n_names)
    for i in 1:n_names
        len = read(io, Int64)
        names[i] = String([read(io, UInt8) for _ in 1:len])
    end
    pool_len = read(io, Int64)
    pool = Vector{UInt64}(undef, pool_len); read!(io, pool)
    n_kmers = length(kct.kmer.seqs)
    ids = Vector{UInt16}(undef, n_kmers); read!(io, ids)
    return KCT(kct.kmer, kct.counts, BiotypLayer(ids, pool, names))
end

# V1.4 (NeoKCT) → KCT{K,Ab,CountsLayer,Nothing,UInt64,UInt32}
function _load_kct(io::IO, ::Val{1.4})
    K = Int(read(io, Int64))
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    n_kmers = read(io, Int64)
    n_flat_cids = read(io, Int64)
    words_len = read(io, Int64)
    bitmap_len = read(io, Int64)
    n_samples = read(io, Int64)
    W = _WORD_TYPES[read(io, Int64)]
    cp_interval = Int(read(io, Int64))
    n_cp = read(io, Int64)
    n_rci = read(io, Int64)
    checkpoints = Vector{UInt64}(undef, n_cp); read!(io, checkpoints)
    deltas = Vector{UInt32}(undef, n_kmers); read!(io, deltas)
    rci_i64 = Vector{Int64}(undef, n_rci); read!(io, rci_i64)
    n_cids = Vector{UInt16}(undef, n_kmers); read!(io, n_cids)
    flat_cids = Vector{UInt32}(undef, n_flat_cids); read!(io, flat_cids)
    words = Vector{W}(undef, words_len); read!(io, words)
    bitmap = BitVector(undef, bitmap_len)
    for i in eachindex(bitmap.chunks); bitmap.chunks[i] = read(io, UInt64); end
    seqs = DeltaArray(checkpoints, deltas, UInt64.(rci_i64), cp_interval)
    kl = KmerLayer{K, Ab, UInt64, UInt32}(seqs, Ref(20) => fill(0:-1, 1 << max(0, (K - 4) * bits_per_symbol(Ab()))))
    cl = CountsLayer(n_cids, flat_cids, PackedArray{UInt32, W}(words, bitmap), Ref(n_samples))
    kct = KCT(kl, cl)
    compute_index!(kct)
    return kct
end

# V1.3 (NeoKCT, flat seqs) → KCT{K,Ab,CountsLayer,Nothing,UInt64,UInt32}
function _load_kct(io::IO, ::Val{1.3})
    K = Int(read(io, Int64))
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    n_kmers = read(io, Int64)
    n_flat_cids = read(io, Int64)
    words_len = read(io, Int64)
    bitmap_len = read(io, Int64)
    n_samples = read(io, Int64)
    W = _WORD_TYPES[read(io, Int64)]
    raw_seqs = Vector{UInt64}(undef, n_kmers); read!(io, raw_seqs)
    n_cids = Vector{UInt16}(undef, n_kmers); read!(io, n_cids)
    flat_cids = Vector{UInt32}(undef, n_flat_cids); read!(io, flat_cids)
    words = Vector{W}(undef, words_len); read!(io, words)
    bitmap = BitVector(undef, bitmap_len)
    for i in eachindex(bitmap.chunks); bitmap.chunks[i] = read(io, UInt64); end
    seqs = DeltaArray(raw_seqs, DEFAULT_CHECKPOINT_INTERVAL)
    kl = KmerLayer{K, Ab, UInt64, UInt32}(seqs, Ref(20) => fill(0:-1, 1 << max(0, (K - 4) * bits_per_symbol(Ab()))))
    cl = CountsLayer(n_cids, flat_cids, PackedArray{UInt32, W}(words, bitmap), Ref(n_samples))
    kct = KCT(kl, cl)
    compute_index!(kct)
    return kct
end

# V1.2 (NeoKCT, flat seqs + CSR offsets) → KCT{K,Ab,CountsLayer,Nothing,UInt64,UInt32}
function _load_kct(io::IO, ::Val{1.2})
    K = Int(read(io, Int64))
    Ab_name_len = read(io, Int64)
    Ab = eval(Symbol(String([read(io, UInt8) for _ in 1:Ab_name_len])))
    n_kmers = read(io, Int64)
    n_flat_cids = read(io, Int64)
    words_len = read(io, Int64)
    bitmap_len = read(io, Int64)
    n_samples = read(io, Int64)
    W = _WORD_TYPES[read(io, Int64)]
    raw_seqs = Vector{UInt64}(undef, n_kmers); read!(io, raw_seqs)
    offsets_u32 = Vector{UInt32}(undef, n_kmers + 1); read!(io, offsets_u32)
    n_cids = UInt16.(diff(offsets_u32))
    flat_cids = Vector{UInt32}(undef, n_flat_cids); read!(io, flat_cids)
    words = Vector{W}(undef, words_len); read!(io, words)
    bitmap = BitVector(undef, bitmap_len)
    for i in eachindex(bitmap.chunks); bitmap.chunks[i] = read(io, UInt64); end
    seqs = DeltaArray(raw_seqs, DEFAULT_CHECKPOINT_INTERVAL)
    kl = KmerLayer{K, Ab, UInt64, UInt32}(seqs, Ref(20) => fill(0:-1, 1 << max(0, (K - 4) * bits_per_symbol(Ab()))))
    cl = CountsLayer(n_cids, flat_cids, PackedArray{UInt32, W}(words, bitmap), Ref(n_samples))
    kct = KCT(kl, cl)
    compute_index!(kct)
    return kct
end

## GENOMIC INDEX RETROCOMPAT ##

function load_gidx(path::String)
    open(path, "r") do io
        version = read(io, Float64)
        return _load_gidx(io, Val(version))
    end
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
    kl = KmerLayer{K, Ab, C_type, D_type}(seqs, Ref(20) => fill(0:-1, 1 << max(0, (K - 4) * bits_per_symbol(Ab()))))
    kct = KCT(kl, BiotypLayer(biotype_ids, pool, biotype_names))
    compute_index!(kct)
    return kct
end
