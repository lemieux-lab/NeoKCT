# Default written format for a V3.0 (CountsLayer) table. A SparseCountsLayer table is
# written as V4.0 instead (see write_kct). Every V1.2..4.0 reader stays available.
const KCT_VERSION = 3.0

const _WORD_TYPES = Dict{Int64, DataType}(
    1 => UInt8, 2 => UInt16, 4 => UInt32, 8 => UInt64, 16 => UInt128
)

## PUBLIC API ##

"""
    write_kct(kct, path)

Serialize any `KCT` to `path` in the versioned binary format. A
`SparseCountsLayer` table is written as V4.0, everything else as V3.0
(`KCT_VERSION`). The leading `Float64` version tag lets `load_kct` pick the
right reader. See the format comment below for the on-disk layout.
"""
function write_kct(kct::KCT{K, Ab, Counts}, path::String) where {K, Ab, Counts}
    ver = Counts === SparseCountsLayer ? 4.0 : KCT_VERSION
    open(path, "w") do io
        write(io, ver)
        _write_kct(io, kct, Val(ver))
    end
end

"""
    load_kct(path) -> KCT

Read a `.kct` file written by `write_kct`, dispatching on its version tag. V3.0
and V4.0 files load directly. V1.2, V1.3, V1.4 (NeoKCT) and V2.0 (RichKCT) are
read through the retrocompat loaders below and come back as a current `KCT`.
The prefix search index is rebuilt on load.
"""
function load_kct(path::String)
    open(path, "r") do io
        version = read(io, Float64)
        return _load_kct(io, Val(version))
    end
end

"""
    get_version(path) -> Float64

Read just the leading version tag of a `.kct` file without loading the table.
"""
function get_version(path::String)
    open(path, "r") do io
        return read(io, Float64)
    end
end

## V3.0 / V4.0 FORMAT ##
# Header: [Float64 version][Int64 K][Int64 Ab_name_len][UInt8... Ab_name][UInt8 layers_mask][Int64 n_kmers]
# KmerLayer: [Int64 sizeof(C)][Int64 sizeof(D)][Int64 cp_interval][Int64 n_cp][Int64 n_rci][C... cps][D... deltas][Int64... rci]
# CountsLayer       (V3.0, mask bit 0): [Int64 n_samples][Int64 n_flat_cids][Int64 words_len][Int64 bitmap_len][Int64 sizeof(W)][UInt16... n_cids][UInt32... flat_cids][W... words][UInt64... bitmap.chunks]
# SparseCountsLayer (V4.0, mask bit 0): [Int64 n_samples][Int64 n_rows][Int64 n_pairs][UInt32... row_id (n_kmers)][UInt64... row_offsets (n_rows+1)][UInt32... row_samples (n_pairs)][UInt32... row_counts (n_pairs)]
# BiotypLayer (mask bit 1): [Int64 n_names]([Int64 len][UInt8... name]...)[Int64 pool_len][UInt64... pool][UInt16... ids]

_layers_mask(::Nothing, ::Nothing) = UInt8(0)
_layers_mask(::CountsLayer, ::Nothing) = UInt8(1)
_layers_mask(::SparseCountsLayer, ::Nothing) = UInt8(1)
_layers_mask(::Nothing, ::BiotypLayer) = UInt8(2)
_layers_mask(::CountsLayer, ::BiotypLayer) = UInt8(3)
_layers_mask(::SparseCountsLayer, ::BiotypLayer) = UInt8(3)

# Shared prefix (version + header + KmerLayer block) for V3.0 and V4.0.
# `string(Ab)` keeps type parameters ("DNAAlphabet{2}"), where `.name.singletonname` drops
# them; `_read_header_and_kmers` parses it back with Meta.parse. For plain `AAAlphabet` both
# forms are the same string, so V3.0 files stay byte-identical.
function _write_header_and_kmers(io::IO, kct::KCT{K, Ab, Counts, Biotype, C, D}) where {K, Ab<:Alphabet, Counts, Biotype, C<:Unsigned, D<:Unsigned}
    Ab_name = string(Ab)
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
end

function _write_kct(io::IO, kct::KCT, ::Val{3.0})
    _write_header_and_kmers(io, kct)
    _write_counts(io, kct.counts)
    _write_biotype(io, kct.biotype)
end

function _write_kct(io::IO, kct::KCT, ::Val{4.0})
    _write_header_and_kmers(io, kct)
    _write_sparse_counts(io, kct.counts)
    _write_biotype(io, kct.biotype)
end

function _write_sparse_counts(io::IO, scl::SparseCountsLayer)
    write(io, Int64(scl.n_samples.x))
    write(io, Int64(length(scl.row_offsets) - 1))  # n_rows
    write(io, Int64(length(scl.row_samples)))  # n_pairs
    write(io, scl.row_id)
    write(io, scl.row_offsets)
    write(io, scl.row_samples)
    write(io, scl.row_counts)
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

# Shared prefix reader for V3.0 and V4.0: returns (layers_mask, n_kmers, KmerLayer).
function _read_header_and_kmers(io::IO)
    K = Int(read(io, Int64))
    Ab_name_len = read(io, Int64)
    Ab = eval(Meta.parse(String([read(io, UInt8) for _ in 1:Ab_name_len])))  # parses "DNAAlphabet{2}" too
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
    kl = KmerLayer{K, Ab, C_type, D_type}(seqs, _empty_kmer_idx(K, Ab(), C_type))
    return layers_mask, n_kmers, kl
end

function _load_kct(io::IO, ::Val{3.0})
    layers_mask, n_kmers, kl = _read_header_and_kmers(io)
    cl = (layers_mask & UInt8(1)) != 0 ? _read_counts(io, n_kmers) : nothing
    bl = (layers_mask & UInt8(2)) != 0 ? _read_biotype(io, n_kmers) : nothing
    kct = _build_kct(kl, cl, bl)
    compute_index!(kct)
    return kct
end

function _load_kct(io::IO, ::Val{4.0})
    layers_mask, n_kmers, kl = _read_header_and_kmers(io)
    cl = (layers_mask & UInt8(1)) != 0 ? _read_sparse_counts(io, n_kmers) : nothing
    bl = (layers_mask & UInt8(2)) != 0 ? _read_biotype(io, n_kmers) : nothing
    kct = _build_kct(kl, cl, bl)
    compute_index!(kct)
    return kct
end

function _read_sparse_counts(io::IO, n_kmers::Int64)
    n_samples = read(io, Int64)
    n_rows = read(io, Int64)
    n_pairs = read(io, Int64)
    row_id = Vector{UInt32}(undef, n_kmers); read!(io, row_id)
    row_offsets = Vector{UInt64}(undef, n_rows + 1); read!(io, row_offsets)
    row_samples = Vector{UInt32}(undef, n_pairs); read!(io, row_samples)
    row_counts = Vector{UInt32}(undef, n_pairs); read!(io, row_counts)
    return SparseCountsLayer(row_id, row_offsets, row_samples, row_counts, Ref(n_samples))
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
_build_kct(kl, scl::SparseCountsLayer, ::Nothing) = KCT(kl, scl)
_build_kct(kl, scl::SparseCountsLayer, bl::BiotypLayer) = KCT(kl, scl, bl)

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
    kl = KmerLayer{K, Ab, UInt64, UInt32}(seqs, _empty_kmer_idx(K, Ab(), UInt64))
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
    kl = KmerLayer{K, Ab, UInt64, UInt32}(seqs, _empty_kmer_idx(K, Ab(), UInt64))
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
    kl = KmerLayer{K, Ab, UInt64, UInt32}(seqs, _empty_kmer_idx(K, Ab(), UInt64))
    cl = CountsLayer(n_cids, flat_cids, PackedArray{UInt32, W}(words, bitmap), Ref(n_samples))
    kct = KCT(kl, cl)
    compute_index!(kct)
    return kct
end

## GENOMIC INDEX RETROCOMPAT ##

"""
    load_gidx(path) -> KCT

Load a legacy `.gidx` genomic-index file (V1.0) and return it as a
`KCT{K, Ab, Nothing, BiotypLayer}`, the same type `build_genomic_index` now
produces. New code should write these with `write_kct` instead.
"""
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
    kl = KmerLayer{K, Ab, C_type, D_type}(seqs, _empty_kmer_idx(K, Ab(), C_type))
    kct = KCT(kl, BiotypLayer(biotype_ids, pool, biotype_names))
    compute_index!(kct)
    return kct
end
