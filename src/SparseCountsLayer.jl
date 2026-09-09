## Sparse Counts Layer (.kct V4.0) ##

"""
    SparseCountsLayer <: AbstractLayer

Inline, block-packed store of per-k-mer sparse count vectors (`.kct` V4.0).

Each k-mer's row is the list of `(sample, count)` pairs for the samples where the k-mer is
present, ascending by sample. Absent samples are implicit and never stored. Rows are packed
inline, in k-mer order, with **no cross-k-mer deduplication**: measurement on real RNA-seq
showed rows longer than four pairs (which hold ~96% of all pairs) are effectively unique, so
a dedup pool plus a per-k-mer `row_id` pointer costs about what it saves and needs a
multi-TB build-time dictionary. Dropping it makes the assembly a single streaming transcode.

K-mers are grouped into fixed blocks of `block_size`. Each block is packed on its own with a
frame-of-reference scheme: one bit width for the row lengths, one for the within-row sample
gaps, one for the counts, plus a short exception list for counts that do not fit 16 bits.
`block_ptr` gives the byte offset of each block in `blob`, so reading k-mer `i` decodes just
its block.

Unlike `CountsLayer` (V3.0: `flat_cids` + `n_cids` + a per-k-mer `PackedArray` word chain
with materialised absent-sample zeros, whose per-k-mer index grew superlinearly with the
sample count), nothing here grows with the sample count except the actual gap and count
widths. Built only by `build_kct_streaming`; a `SparseCountsLayer` KCT is read-only
(`push!` / `collapse!` / `repack` / `sort!` all error).

# Fields
- `block_size`: k-mers per block
- `block_ptr`: 0-based byte offset of each block in `blob`, `length` is the block count
- `blob`: the packed blocks, concatenated
- `n_kmers`: total k-mers (row data is inline, so this is not derivable from a pool length)
- `n_samples`: total samples represented (sizes the zero-padded read)
"""
mutable struct SparseCountsLayer <: AbstractLayer
    block_size::Int32
    block_ptr::Vector{UInt64}
    blob::Vector{UInt8}
    n_kmers::Base.RefValue{Int64}
    n_samples::Base.RefValue{Int64}
end

SparseCountsLayer(block_size::Integer, block_ptr::Vector{UInt64}, blob::Vector{UInt8},
                  n_kmers::Integer, n_samples::Integer) =
    SparseCountsLayer(Int32(block_size), block_ptr, blob, Ref(Int64(n_kmers)), Ref(Int64(n_samples)))

Base.length(scl::SparseCountsLayer) = Int(scl.n_kmers.x)

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
function _decode_block(scl::SparseCountsLayer, b::Int)
    blob = scl.blob
    base = Int(scl.block_ptr[b + 1])
    @inbounds nk    = Int(blob[base + 1])
    @inbounds lbits = Int(blob[base + 2])
    @inbounds gbits = Int(blob[base + 3])
    @inbounds cbits = Int(blob[base + 4])
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

# One k-mer's (samples, counts). Decodes its whole block; fine for the scattered lookups the
# walk does. A full-table pass should iterate blocks via `_each_row`.
function _decode_row(scl::SparseCountsLayer, i::Integer)
    B = Int(scl.block_size)
    b = (Int(i) - 1) ÷ B
    li = (Int(i) - 1) % B
    return _decode_block(scl, b)[li + 1]
end

# f(i, samples, counts) for every k-mer in order, decoding each block once.
function _each_row(f, scl::SparseCountsLayer)
    B = Int(scl.block_size)
    n = length(scl)
    i = 0
    for b in 0:(length(scl.block_ptr) - 1)
        for (s, c) in _decode_block(scl, b)
            i += 1
            i > n && return
            f(i, s, c)
        end
    end
end

# Full length-`n_samples` count vector for k-mer `i`: scatter its row into a zero vector.
# Same contract as the CountsLayer method in KCTLayers.jl, always a `Vector{UInt32}`.
function assemble_count_vector(scl::SparseCountsLayer, i::Integer)
    s, c = _decode_row(scl, i)
    v = zeros(UInt32, Int(scl.n_samples.x))
    @inbounds for t in eachindex(s)
        v[s[t]] = c[t]
    end
    return v
end

Base.getindex(scl::SparseCountsLayer, i::Integer) = assemble_count_vector(scl, i)
