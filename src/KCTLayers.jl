using Kmers
using BioSequences
using BioSymbols
using ProgressMeter
using Base.Threads
using NArrays
using Dates
using BitIntegers
using Mmap

abstract type AbstractLayer end

## K-mer Layer ##

const KmerIdxEntry{C} = Tuple{C, UnitRange{Int64}}

"""
    KmerLayer{K, Ab, C, D} <: AbstractLayer

Sorted, delta-compressed store of all unique k-mers in the table.

K-mers are stored as their raw `UInt64` bit-encoding in a `DeltaArray`, which
exploits the sorted order to keep only small deltas between consecutive entries.
A prefix-based binary search index (`idx`) partitions the space by the leading
`prefix_size` symbols so lookups skip directly to the right region.

# Type Parameters
- `K`: k-mer length in amino-acid symbols
- `Ab`: alphabet type (e.g. `AAAlphabet`)
- `C`: checkpoint word type for the `DeltaArray` (default `UInt64`)
- `D`: delta word type for the `DeltaArray` (default `UInt32`)

# Fields
- `seqs`: delta-compressed sorted k-mer bit-encodings
- `idx`: prefix-length reference paired with a per-prefix `(value, range)` lookup table
  `value` is the decoded k-mer value at `range.start`, letting findfirst seed searchfirst
  directly instead of rewinding to the nearest DeltaArray checkpoint (see the 5-arg
  `searchfirst` in NArrays)
"""
struct KmerLayer{K, Ab <: Alphabet, C <: Unsigned, D <: Unsigned} <: AbstractLayer
    seqs::DeltaArray{C, D}
    idx::Pair{Base.RefValue{Int64}, Vector{KmerIdxEntry{C}}}
end

idx_prefix_size(klayer::KmerLayer) = klayer.idx[1].x

function Base.getindex(klayer::KmerLayer, i::UnitRange)
    return Tuple(klayer[j] for j in i)
end

function Base.getindex(klayer::KmerLayer, i::AbstractVector)
    return Tuple(klayer[j] for j in i)
end

const DEFAULT_CHECKPOINT_INTERVAL = 256

Base.length(kl::KmerLayer) = length(kl.seqs)

function Base.getindex(kl::KmerLayer{K, Ab}, i::Integer) where {K, Ab<:Alphabet}
    return Kmer{Ab, K, 1}(Kmers.unsafe, (kl.seqs[i],))
end

# Prefix bits the index does NOT partition on; trades index RAM against per-lookup scan
# length (see ARCHITECTURE.md). 5 keeps the AA index sub-GB with a negligible scan cost.
const DEFAULT_IDX_PREFIX_SIZE = 5

# Empty prefix index sized for a K-mer table over alphabet Ab, with checkpoint word type C.
_empty_kmer_idx(K::Integer, ab::Alphabet, ::Type{C};
                prefix_size::Integer=DEFAULT_IDX_PREFIX_SIZE) where {C<:Unsigned} =
    Ref(Int64(prefix_size)) => fill((zero(C), 0:-1),
                                    1 << max(0, (K - prefix_size) * bits_per_symbol(ab)))

KmerLayer{K, Ab}(; checkpoint_size::Type{C}=UInt64,
                   delta_size::Type{D}=UInt32) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KmerLayer{K, Ab, C, D}(DeltaArray{checkpoint_size, delta_size}(DEFAULT_CHECKPOINT_INTERVAL),
                            _empty_kmer_idx(K, Ab(), C))

function Base.findfirst(kl::KmerLayer{K, Ab}, key::Kmer{Ab, K}) where {K, Ab<:Alphabet}
    key_bits = key.data[1]
    idx_key = (key_bits >> (idx_prefix_size(kl) * bits_per_symbol(Ab()))) + 1
    lo_val, r = kl.idx[2][idx_key]
    isempty(r) && return 0
    return searchfirst(kl.seqs, key_bits, r.start, r.stop, lo_val)
end

Base.findfirst(kl::KmerLayer{K, Ab}, key::UInt64) where {K, Ab<:Alphabet} =
    findfirst(kl, Kmer{Ab, K, 1}(Kmers.unsafe, (key,)))

"""
    compute_index!(kl; prefix_size=DEFAULT_IDX_PREFIX_SIZE)

Rebuild the prefix search index of a `KmerLayer` (or a `KCT`, which forwards to
its `kmer` layer). Walks the sorted k-mer array once and records, for each
distinct leading-`prefix_size`-symbol group, the `(decoded value at group start,
index range)` pair that `findfirst` uses to seed `searchfirst` directly. Call it
after any operation that changes the k-mer set. The layer's index vector must
already be sized for `prefix_size` (see `_empty_kmer_idx`).
"""
function compute_index!(kl::KmerLayer{K, Ab, C};
                        prefix_size::Int64=DEFAULT_IDX_PREFIX_SIZE) where {K, Ab<:Alphabet, C<:Unsigned}
    n = length(kl.seqs)
    want = 1 << max(0, (K - prefix_size) * bits_per_symbol(Ab()))
    @assert length(kl.idx[2]) == want "prefix index sized for a different prefix_size " *
        "($(length(kl.idx[2])) buckets, expected $want): rebuild the KmerLayer via _empty_kmer_idx(...; prefix_size)"
    kl.idx[1].x = prefix_size
    n == 0 && return
    prefix_shift = prefix_size * bits_per_symbol(Ab())

    start = 1
    start_val = zero(C)
    last_key = zero(C)
    seeded = false  # start_val isn't real until the first (i=1) iteration runs

    # Record where each distinct prefix group starts and its decoded start value, so
    # findfirst can seed searchfirst directly. `seeded` guards the first group (i=1): sealing
    # it with the zero(C) placeholder would misfire if a real k-mer decodes to exactly 0
    # (see ARCHITECTURE.md).
    @showprogress "Computing Binary Search Index..." for (i, val) in enumerate(kl.seqs)
        key = val >> prefix_shift
        if !seeded || key > last_key
            seeded && (kl.idx[2][last_key+1] = (start_val, start:i))
            start = i
            start_val = val
            last_key = key
            seeded = true
        end
    end
    kl.idx[2][last_key+1] = (start_val, start:n)
end

## Counts Layer ##

"""
    CountsLayer <: AbstractLayer

Inline, block-packed store of per-k-mer sparse count vectors (`.kct` V4.0). Each k-mer's row
is its `(sample, count)` pairs, ascending by sample, absent samples implicit. Rows are packed
inline in k-mer order with no cross-k-mer deduplication (see ARCHITECTURE.md for why).
K-mers are grouped into fixed blocks of `block_size`, each packed with a frame-of-reference
scheme (one bit width per field, short exception list for large counts); `block_ptr` gives
each block's byte offset in `blob` for direct access. Built only by `build_kct_streaming`;
read-only (`push!` / `collapse!` / `repack` / `sort!` all error).

# Fields
- `block_size`: k-mers per block
- `block_ptr`: 0-based byte offset of each block in `blob`, `length` is the block count
- `blob`: the packed blocks, concatenated
- `n_kmers`: total k-mers (row data is inline, so this is not derivable from a pool length)
- `n_samples`: total samples represented (sizes the zero-padded read)
"""
mutable struct CountsLayer <: AbstractLayer
    block_size::Int32
    block_ptr::Vector{UInt64}  # k-mer i's block at block_ptr[i/block_size]
    blob::Vector{UInt8}
    n_kmers::Base.RefValue{Int64}
    n_samples::Base.RefValue{Int64}
end

CountsLayer(block_size::Integer, block_ptr::Vector{UInt64}, blob::Vector{UInt8},
            n_kmers::Integer, n_samples::Integer) =
    CountsLayer(Int32(block_size), block_ptr, blob, Ref(Int64(n_kmers)), Ref(Int64(n_samples)))

Base.length(cl::CountsLayer) = Int(cl.n_kmers.x)

## bit + integer stream codec ##

# MSB-first bit packer that appends whole bytes to an existing Vector{UInt8}. `bits` per
# call is small here (<= 16), so the accumulator never holds more than ~23 live bits.
mutable struct _BitWriter
    buf::Vector{UInt8}
    acc::UInt64
    nb::Int
end
_BitWriter(buf::Vector{UInt8}) = _BitWriter(buf, UInt64(0), 0)

@inline function _put!(w::_BitWriter, v::UInt64, bits::Int)
    bits == 0 && return
    mask = bits == 64 ? typemax(UInt64) : (UInt64(1) << bits) - UInt64(1)
    w.acc = (w.acc << bits) | (v & mask)
    w.nb += bits
    while w.nb >= 8
        w.nb -= 8
        push!(w.buf, UInt8((w.acc >> w.nb) & 0xff))
    end
end

@inline function _flush!(w::_BitWriter)
    if w.nb > 0
        push!(w.buf, UInt8((w.acc << (8 - w.nb)) & 0xff))
        w.nb = 0
        w.acc = UInt64(0)
    end
end

mutable struct _BitReader
    buf::Vector{UInt8}
    byte::Int          # 1-based index of the next unread byte
    acc::UInt64
    nb::Int
end
_BitReader(buf::Vector{UInt8}, byte::Int) = _BitReader(buf, byte, UInt64(0), 0)

@inline function _get!(r::_BitReader, bits::Int)::UInt64
    bits == 0 && return UInt64(0)
    while r.nb < bits
        @inbounds r.acc = (r.acc << 8) | UInt64(r.buf[r.byte])
        r.byte += 1
        r.nb += 8
    end
    r.nb -= bits
    mask = bits == 64 ? typemax(UInt64) : (UInt64(1) << bits) - UInt64(1)
    return (r.acc >> r.nb) & mask
end

@inline function _byte_align!(r::_BitReader)
    r.nb = 0
    r.acc = UInt64(0)
end

@inline _bitlen(x::Integer) = x <= 0 ? 0 : Base.top_set_bit(unsigned(x))

@inline function _push_u32le!(buf::Vector{UInt8}, v::UInt32)
    push!(buf, UInt8(v & 0xff), UInt8((v >> 8) & 0xff), UInt8((v >> 16) & 0xff), UInt8((v >> 24) & 0xff))
end
@inline function _get_u32le(buf::Vector{UInt8}, byte::Int)
    @inbounds v = UInt32(buf[byte]) | (UInt32(buf[byte + 1]) << 8) |
                  (UInt32(buf[byte + 2]) << 16) | (UInt32(buf[byte + 3]) << 24)
    return v, byte + 4
end

function _put_varint!(buf::Vector{UInt8}, x::UInt64)
    while true
        b = UInt8(x & 0x7f)
        x >>= 7
        if x != 0
            push!(buf, b | 0x80)
        else
            push!(buf, b)
            return
        end
    end
end
function _get_varint(buf::Vector{UInt8}, byte::Int)
    x = UInt64(0); shift = 0
    while true
        @inbounds b = buf[byte]
        byte += 1
        x |= UInt64(b & 0x7f) << shift
        (b & 0x80) == 0 && return x, byte
        shift += 7
    end
end

## block packing ##

# Block layout in `blob` (starts at block_ptr[b]+1, 1-based):
#   UInt8 nk | UInt8 lbits | UInt8 gbits | UInt8 cbits | varint n_exc
#   bitstream (MSB-first): nk lengths (lbits) | T gaps (gbits) | T counts (cbits), T = sum(lengths)
#   byte-aligned, then n_exc records: varint local_pair_index | UInt32LE value
# lbits==0 means every length is 1; cbits==0 means every count is 1. A packed count of
# 0xffff with n_exc>0 is a sentinel resolved from the exception list.
const _V4_COUNT_INLINE_MAX = UInt32(0xffff)

# Rows come in flat CSR: row j (1..nk) is `s[off[j]:off[j+1]-1]` / `c[...]`, with off[1] the
# base into s/c. Fed straight from the assembly's merge buffers so nothing allocates a
# Vector per k-mer.
function _pack_block!(blob::Vector{UInt8}, s::Vector{UInt32}, c::Vector{UInt32},
                      off::Vector{Int}, nk::Int)
    @assert 1 <= nk <= 255
    base = off[1] - 1                # pair p of the block is s[base + p], p in 1:T
    T = off[nk + 1] - off[1]

    maxlen = 0
    @inbounds for j in 1:nk
        L = off[j + 1] - off[j]
        L > maxlen && (maxlen = L)
    end
    lbits = maxlen <= 1 ? 0 : _bitlen(maxlen)

    gaps = Vector{UInt32}(undef, T)
    maxgap = UInt32(1); maxc = UInt32(0)
    @inbounds for j in 1:nk
        prev = UInt32(0)
        for p in off[j]:(off[j + 1] - 1)
            g = s[p] - prev; prev = s[p]
            gaps[p - base] = g
            g > maxgap && (maxgap = g)
            cc = c[p]
            cc > maxc && (maxc = cc)
        end
    end
    gbits = max(1, _bitlen(maxgap))

    exc = Tuple{Int, UInt32}[]
    if maxc < _V4_COUNT_INLINE_MAX
        cbits = maxc <= 1 ? 0 : _bitlen(maxc)
    else
        cbits = 16
        @inbounds for p in 1:T
            c[base + p] >= _V4_COUNT_INLINE_MAX && push!(exc, (p - 1, c[base + p]))
        end
    end

    push!(blob, UInt8(nk), UInt8(lbits), UInt8(gbits), UInt8(cbits))
    _put_varint!(blob, UInt64(length(exc)))

    w = _BitWriter(blob)
    if lbits > 0
        @inbounds for j in 1:nk
            _put!(w, UInt64(off[j + 1] - off[j]), lbits)
        end
    end
    @inbounds for p in 1:T
        _put!(w, UInt64(gaps[p]), gbits)
    end
    if cbits > 0
        @inbounds for p in 1:T
            cc = c[base + p]
            _put!(w, UInt64(cc >= _V4_COUNT_INLINE_MAX ? _V4_COUNT_INLINE_MAX : cc), cbits)
        end
    end
    _flush!(w)

    for (idx, val) in exc
        _put_varint!(blob, UInt64(idx))
        _push_u32le!(blob, val)
    end
    return nothing
end

# Decode every row of block `b` (0-based): returns nk `(samples, counts)` pairs.
function _decode_block(cl::CountsLayer, b::Int)
    blob = cl.blob
    base = Int(cl.block_ptr[b + 1])
    @inbounds nk = Int(blob[base + 1])
    @inbounds lbits = Int(blob[base + 2])  # bit width shared by every row length in the block
    @inbounds gbits = Int(blob[base + 3])  # bit width shared by every sample-gap in the block
    @inbounds cbits = Int(blob[base + 4])  # bit width shared by every count in the block (capped at 16, exc otherwise)
    n_exc, byte = _get_varint(blob, base + 5)

    r = _BitReader(blob, byte)
    lens = Vector{Int}(undef, nk)
    if lbits == 0
        fill!(lens, 1)
    else
        @inbounds for j in 1:nk
            lens[j] = Int(_get!(r, lbits))
        end
    end
    T = 0
    @inbounds for j in 1:nk
        T += lens[j]
    end

    gaps = Vector{UInt32}(undef, T)
    @inbounds for i in 1:T
        gaps[i] = UInt32(_get!(r, gbits))
    end
    cnts = Vector{UInt32}(undef, T)
    if cbits == 0
        fill!(cnts, UInt32(1))
    else
        @inbounds for i in 1:T
            cnts[i] = UInt32(_get!(r, cbits))
        end
    end
    if n_exc > 0
        _byte_align!(r)
        eb = r.byte
        for _ in 1:n_exc
            idx, eb = _get_varint(blob, eb)
            val, eb = _get_u32le(blob, eb)
            @inbounds cnts[Int(idx) + 1] = val
        end
    end

    out = Vector{Tuple{Vector{UInt32}, Vector{UInt32}}}(undef, nk)
    off = 0
    @inbounds for j in 1:nk
        L = lens[j]
        s = Vector{UInt32}(undef, L)
        acc = UInt32(0)
        for t in 1:L
            acc += gaps[off + t]
            s[t] = acc
        end
        out[j] = (s, cnts[off + 1:off + L])
        off += L
    end
    return out
end

# One k-mer's (samples, counts). Decodes its whole block, fine for the scattered lookups the
# walk does. A full-table pass should iterate blocks via `_each_row`.
function _decode_row(cl::CountsLayer, i::Integer)
    B = Int(cl.block_size)
    b = (Int(i) - 1) ÷ B
    li = (Int(i) - 1) % B
    return _decode_block(cl, b)[li + 1]
end

# f(i, samples, counts) for every k-mer in order, decoding each block once.
function _each_row(f, cl::CountsLayer)
    B = Int(cl.block_size)
    n = length(cl)
    i = 0
    for b in 0:(length(cl.block_ptr) - 1)
        for (s, c) in _decode_block(cl, b)
            i += 1
            i > n && return
            f(i, s, c)
        end
    end
end

# Full length-`n_samples` count vector for k-mer `i`: scatter its row into a zero vector.
function assemble_count_vector(cl::CountsLayer, i::Integer)
    s, c = _decode_row(cl, i)
    v = zeros(UInt32, Int(cl.n_samples.x))
    @inbounds for t in eachindex(s)
        v[s[t]] = c[t]
    end
    return v
end

Base.getindex(cl::CountsLayer, i::Integer) = assemble_count_vector(cl, i)

## Biotype Layer ##

const INTERGENIC_MASK = UInt64(1)  # bit 0 reserved for intergenic

"""
    BiotypLayer <: AbstractLayer

Per-k-mer biotype annotation stored as interned bitmasks.

Each bit position corresponds to one biotype in `biotype_names`: bit `b-1` set
means the k-mer originates from `biotype_names[b]`. Bitmasks are deduplicated
into a `pool` so `ids` stores only a `UInt16` index per k-mer. `pool[1]` is
always `INTERGENIC_MASK` (bit 0 only), and k-mers with no annotation default to it.

# Fields
- `ids`: one pool index per k-mer (1-based into `pool`)
- `pool`: deduplicated bitmasks, `pool[1] == INTERGENIC_MASK`
- `biotype_names`: human-readable name per bit position (`biotype_names[b]` ↔ bit `b-1`)
"""
struct BiotypLayer <: AbstractLayer
    ids::Vector{UInt16}  # one per k-mer, index into pool (1-based)
    pool::Vector{UInt64}  # deduplicated biotype bitmasks, pool[1] always == INTERGENIC_MASK
    biotype_names::Vector{String}  # biotype_names[b] = name for bit b-1
end

# All-intergenic layer. biotype_names[1] must be "intergenic".
BiotypLayer(n::Int, biotype_names::Vector{String}) =
    BiotypLayer(fill(UInt16(1), n), UInt64[INTERGENIC_MASK], biotype_names)

Base.length(bl::BiotypLayer) = length(bl.ids)

biotype_mask(l::BiotypLayer, i::Int) = l.pool[l.ids[i]]

function biotype_names_for(l::BiotypLayer, i::Int)
    mask = biotype_mask(l, i)
    return [l.biotype_names[b] for b in eachindex(l.biotype_names) if (mask >> (b - 1)) & 1 == 1]
end

function has_biotype(l::BiotypLayer, i::Int, name::String)
    b = findfirst(==(name), l.biotype_names)
    isnothing(b) && throw(ArgumentError("Unknown biotype: $name"))
    return (biotype_mask(l, i) >> (b - 1)) & 1 == 1
end

# Insert mask into pool if absent, then return its 1-based UInt16 index.
function _intern_mask!(pool::Vector{UInt64}, index::Dict{UInt64, UInt16}, mask::UInt64)::UInt16
    return get!(index, mask) do
        push!(pool, mask)
        UInt16(length(pool))
    end
end

## KCT Wrapper ##

"""
    KCT{K, Ab, Counts, Biotype, C, D}

K-mer Count Table. Top-level container combining a sorted k-mer store with
optional per-sample count and biotype-annotation layers.

The table is always keyed by a `KmerLayer` storing k-mers in sorted order with
O(log n) prefix-indexed lookup. `Counts` and `Biotype` are either the
corresponding layer type or `Nothing`, letting the compiler dispatch to
specialised `getindex` methods with no runtime overhead.

# Type Parameters
- `K`: k-mer length in amino-acid symbols
- `Ab`: alphabet type (e.g. `AAAlphabet`)
- `Counts`: `CountsLayer` when count data is present, `Nothing` otherwise
- `Biotype`: `BiotypLayer` when biotype data is present, `Nothing` otherwise
- `C`: checkpoint word type forwarded to the underlying `KmerLayer`
- `D`: delta word type forwarded to the underlying `KmerLayer`

# Fields
- `kmer`: sorted k-mer store
- `counts`: per-sample count layer, or `nothing`
- `biotype`: per-k-mer biotype annotation layer, or `nothing`
"""
struct KCT{K, Ab <: Alphabet, Counts <: Union{CountsLayer, Nothing}, Biotype <: Union{BiotypLayer, Nothing}, C <: Unsigned, D <: Unsigned}
    kmer::KmerLayer{K, Ab, C, D}
    counts::Counts
    biotype::Biotype
end

KCT(kl::KmerLayer{K, Ab, C, D}) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, Nothing, Nothing, C, D}(kl, nothing, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, Nothing, BiotypLayer, C, D}(kl, nothing, bl)
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, Nothing, C, D}(kl, cl, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, BiotypLayer, C, D}(kl, cl, bl)

## KCT Methods ##

Base.length(kct::KCT) = length(kct.kmer)
Base.size(kct::KCT) = (length(kct),)
idx_prefix_size(kct::KCT) = idx_prefix_size(kct.kmer)
compute_index!(kct::KCT; kwargs...) = compute_index!(kct.kmer; kwargs...)

Base.findfirst(kct::KCT{K, Ab}, key::Kmer{Ab, K}) where {K, Ab} = findfirst(kct.kmer, key)
Base.findfirst(kct::KCT{K, Ab}, key::UInt64) where {K, Ab} = findfirst(kct.kmer, key)

function Base.getindex(kct::KCT{K, Ab, Nothing, Nothing}, i::Integer) where {K, Ab}
    return kct.kmer[i]
end

function Base.getindex(kct::KCT{K, Ab, Nothing, BiotypLayer}, i::Integer) where {K, Ab}
    return kct.kmer[i] => biotype_mask(kct.biotype, i)
end

function Base.getindex(kct::KCT{K, Ab, CountsLayer, Nothing}, i::Integer) where {K, Ab}
    return kct.kmer[i] => kct.counts[i]
end

function Base.getindex(kct::KCT{K, Ab, CountsLayer, BiotypLayer}, i::Integer) where {K, Ab}
    return kct.kmer[i] => (; counts=kct.counts[i], biotype=biotype_mask(kct.biotype, i))
end

Base.getindex(kct::KCT, i::UnitRange) = Tuple(kct[j] for j in i)
Base.getindex(kct::KCT, i::AbstractVector) = Tuple(kct[j] for j in i)

# A CountsLayer table is assembled whole by build_kct_streaming and never mutated afterwards.
const _COUNTS_READONLY = "CountsLayer KCTs are read-only; build them with build_kct_streaming"
collapse!(::KCT{K, Ab, CountsLayer}) where {K, Ab} = error(_COUNTS_READONLY)
Base.sort!(::KCT{K, Ab, CountsLayer}) where {K, Ab} = error(_COUNTS_READONLY)
Base.push!(::KCT{K, Ab, CountsLayer}, ::Any) where {K, Ab} = error(_COUNTS_READONLY)
repack(::KCT{K, Ab, CountsLayer}, ::Type{<:Unsigned}; kwargs...) where {K, Ab} = error(_COUNTS_READONLY)

## KCT Genomic Index ##

# Build a KCT{K, Ab, Nothing, BiotypLayer} from a sorted k-mer vector and a parallel bitmask vector.
# Replaces the old GenomicIndex{K,Ab} constructor.
function KCT{K, Ab}(sorted_kmers::Vector{UInt64}, bitmasks::Vector{UInt64},
                    biotype_names::Vector{String};
                    checkpoint_size::Type{C}=UInt64,
                    delta_size::Type{D}=UInt32) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned}
    @assert length(sorted_kmers) == length(bitmasks)
    @assert !isempty(biotype_names) && biotype_names[1] == "intergenic"
    pool = UInt64[]
    index_map = Dict{UInt64, UInt16}()
    ids = Vector{UInt16}(undef, length(bitmasks))
    for (i, mask) in enumerate(bitmasks)
        ids[i] = _intern_mask!(pool, index_map, mask)
    end
    kl = KmerLayer{K, Ab, C, D}(DeltaArray{checkpoint_size, delta_size}(sorted_kmers, DEFAULT_CHECKPOINT_INTERVAL),
                                  _empty_kmer_idx(K, Ab(), C))
    kct = KCT(kl, BiotypLayer(ids, pool, biotype_names))
    compute_index!(kct.kmer)
    return kct
end

"""
    add_biotypes(kct, gidx) -> KCT

Left-join the k-mers of `kct` against a genomic index `gidx`
(`KCT{K, Ab, Nothing, BiotypLayer}`, from `build_genomic_index`) and return a new
`KCT` carrying a `BiotypLayer`. The join is an O(n + m) walk over the two sorted
k-mer stores. A k-mer that matches an entry in `gidx` takes that entry's biotype
mask. A k-mer with no match keeps the intergenic mask.

The body only walks `kct.kmer.seqs` and forwards `kct.counts` untouched, so it
works whether or not `kct` carries a `CountsLayer`.
"""
function add_biotypes(kct::KCT{K, Ab, Counts, Nothing, C, D},
                      gidx::KCT{K, Ab, Nothing, BiotypLayer}) where {K, Ab, Counts<:Union{CountsLayer, Nothing}, C<:Unsigned, D<:Unsigned}
    n = length(kct.kmer)
    pool = UInt64[INTERGENIC_MASK]
    index_map = Dict{UInt64, UInt16}(INTERGENIC_MASK => UInt16(1))
    ids = fill(UInt16(1), n)

    kct_iter = iterate(kct.kmer.seqs)
    gidx_iter = iterate(gidx.kmer.seqs)
    i = 1; j = 1

    # O(n + m) sorted merge walk: advance the smaller pointer each step,
    # assign a biotype mask only on exact k-mer match.
    while !isnothing(kct_iter) && !isnothing(gidx_iter)
        kct_val, kct_state = kct_iter
        gidx_val, gidx_state = gidx_iter
        if kct_val == gidx_val
            ids[i] = _intern_mask!(pool, index_map, biotype_mask(gidx.biotype, j))
            kct_iter = iterate(kct.kmer.seqs, kct_state)
            gidx_iter = iterate(gidx.kmer.seqs, gidx_state)
            i += 1; j += 1
        elseif kct_val < gidx_val
            kct_iter = iterate(kct.kmer.seqs, kct_state)
            i += 1
        else
            gidx_iter = iterate(gidx.kmer.seqs, gidx_state)
            j += 1
        end
    end

    n_intergenic = count(==(UInt16(1)), ids)
    printstyled("Biotype assignment done: $n_intergenic / $n intergenic, $(length(pool)) unique bitmasks\n",
                color=:green)
    return KCT(kct.kmer, kct.counts, BiotypLayer(ids, pool, gidx.biotype.biotype_names))
end

include("JelloFish.jl")
include("JellyfishDump.jl")
include("GenomicIndexBuilder.jl")
include("KCTBenchmarker.jl")
include("KCTLoader.jl")

# From-scratch streaming builder: k-way merge of per-sample sorted k-mer streams into a
# CountsLayer table. Needs load_kct / write_kct from KCTLoader.jl above.
include("StreamBuild.jl")
