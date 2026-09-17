# The only version this writes and reads is V4.0 (CountsLayer). V1.2-3.0 (the old
# incremental CountsLayer/RichKCT formats) are gone as of the CountsLayer rename; read one
# with a checkout of the V4.0-last-retrocompat tag.
const KCT_VERSION = 4.0

const _WORD_TYPES = Dict{Int64, DataType}(
    1 => UInt8, 2 => UInt16, 4 => UInt32, 8 => UInt64, 16 => UInt128
)

## PUBLIC API ##

"""
    write_kct(kct, path)

Serialize a `KCT` to `path` in the versioned binary format. The leading
`Float64` version tag lets `load_kct` pick the right reader. See the format
comment below for the on-disk layout.
"""
function write_kct(kct::KCT, path::String)
    open(path, "w") do io
        write(io, KCT_VERSION)
        _write_kct(io, kct, Val(KCT_VERSION))
    end
end

"""
    load_kct(path) -> KCT

Read a `.kct` file written by `write_kct`, dispatching on its version tag. The
prefix search index is rebuilt on load.
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

## V4.0 FORMAT ##
# Header: [Float64 version][Int64 K][Int64 Ab_name_len][UInt8... Ab_name][UInt8 layers_mask][Int64 n_kmers]
# KmerLayer: [Int64 sizeof(C)][Int64 sizeof(D)][Int64 cp_interval][Int64 n_cp][Int64 n_rci][C... cps][D... deltas][Int64... rci]
# CountsLayer (mask bit 0): [Int64 n_samples][Int64 n_kmers][Int32 block_size][Int64 n_blocks][Int64 blob_len][UInt64... block_ptr (n_blocks)][UInt8... blob (blob_len)]
# BiotypLayer (mask bit 1): [Int64 n_names]([Int64 len][UInt8... name]...)[Int64 pool_len][UInt64... pool][UInt16... ids]

_layers_mask(::Nothing, ::Nothing) = UInt8(0)
_layers_mask(::CountsLayer, ::Nothing) = UInt8(1)
_layers_mask(::Nothing, ::BiotypLayer) = UInt8(2)
_layers_mask(::CountsLayer, ::BiotypLayer) = UInt8(3)

# `string(Ab)` keeps type parameters ("DNAAlphabet{2}"), where `.name.singletonname` drops
# them; `_read_header_and_kmers` parses it back with Meta.parse. For plain `AAAlphabet` both
# forms are the same string.
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

function _write_kct(io::IO, kct::KCT, ::Val{4.0})
    _write_header_and_kmers(io, kct)
    _write_counts(io, kct.counts)
    _write_biotype(io, kct.biotype)
end

# Counts-section header: everything up to but not including the blob bytes. Split out so
# _assemble_v4 can write the header and then stream a multi-hundred-GB blob straight from a
# temp file instead of holding it in RAM.
function _write_counts_header(io::IO, n_samples::Integer, n_kmers::Integer, block_size::Integer,
                              n_blocks::Integer, blob_len::Integer, block_ptr::Vector{UInt64})
    write(io, Int64(n_samples))
    write(io, Int64(n_kmers))
    write(io, Int32(block_size))
    write(io, Int64(n_blocks))
    write(io, Int64(blob_len))
    write(io, block_ptr)
end

_write_counts(::IO, ::Nothing) = nothing
function _write_counts(io::IO, cl::CountsLayer)
    _write_counts_header(io, cl.n_samples.x, cl.n_kmers.x, cl.block_size,
                         length(cl.block_ptr), length(cl.blob), cl.block_ptr)
    write(io, cl.blob)
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

# Returns (layers_mask, n_kmers, KmerLayer).
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

function _load_kct(io::IO, ::Val{4.0})
    layers_mask, n_kmers, kl = _read_header_and_kmers(io)
    cl = (layers_mask & UInt8(1)) != 0 ? _read_counts(io, n_kmers) : nothing
    bl = (layers_mask & UInt8(2)) != 0 ? _read_biotype(io, n_kmers) : nothing
    kct = _build_kct(kl, cl, bl)
    compute_index!(kct)
    return kct
end

function _read_counts(io::IO, n_kmers::Int64)
    n_samples = read(io, Int64)
    nk = read(io, Int64)                       # == n_kmers (header), kept for self-consistency
    block_size = read(io, Int32)
    n_blocks = read(io, Int64)
    blob_len = read(io, Int64)
    block_ptr = Vector{UInt64}(undef, n_blocks); read!(io, block_ptr)
    # mmap the blob rather than reading it in: at full cohort it is ~1-2 TB, and a query
    # touches one block. The mapping outlives `io` being closed. Falls back to a plain read
    # if mmap is unavailable (e.g. a filesystem that refuses it).
    blob = if blob_len == 0
        UInt8[]
    else
        try
            b = Mmap.mmap(io, Vector{UInt8}, Int(blob_len))
            seek(io, position(io) + blob_len)
            b
        catch
            b = Vector{UInt8}(undef, blob_len); read!(io, b); b
        end
    end
    return CountsLayer(block_size, block_ptr, blob, nk, n_samples)
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
