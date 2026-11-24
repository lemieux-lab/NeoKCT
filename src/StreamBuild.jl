## Streaming k-way merge builder -> .kct V4.0 (CountsLayer) ##
#
# External merge sort for cohorts too large for push!/collapse! (see ARCHITECTURE.md for the
# full pipeline: wave count -> shard scatter -> fold -> interleaved hierarchical merge ->
# adaptive reshard -> streaming assembly). Alphabet-generic below the counter; only the
# shard shift and the final KmerLayer need `translate` (DNA -> KCT{K,DNAAlphabet{2}}, AA ->
# KCT{K÷3,AAAlphabet}).
const _STREAM_BITS_PER_AA = 5   # bits_per_symbol(AAAlphabet())
const _STREAM_BITS_PER_NT = 2   # bits_per_symbol(DNAAlphabet{2}())

# _log is defined in Progress.jl (included earlier), shared with _Prog/tick!.

# One (k-mer, sample, count) observation, written raw to the wave shard record files.
struct _StreamRec
    kmer::UInt64
    slot::UInt32  # 1-based GLOBAL sample index
    count::UInt32
end
Base.isless(a::_StreamRec, b::_StreamRec) = isless(a.kmer, b.kmer)

# Shard index (0-based) of a k-mer: its top `shard_bits` encoded bits. `mask` = P-1 guards
# against a stray high bit landing the record outside the shard array.
@inline _shard_of(kbits::UInt64, keep_shift::Int, mask::UInt64) = Int((kbits >> keep_shift) & mask)

# One prefix shard of a partial table: distinct k-mers (sorted) each pointing via `row_id`
# into a shard-local pool of sparse rows (`offsets` CSR over `pool_s` / `pool_c`).
struct _PartialShard
    kmers::Vector{UInt64}
    row_id::Vector{UInt32}
    offsets::Vector{UInt64}  # length n_rows + 1, 1-based, offsets[1] == 1
    pool_s::Vector{UInt32}
    pool_c::Vector{UInt32}
end

## Partial-shard disk IO ##

#= On-disk partial shards hold every (sample, count) pair of the cohort unpacked, which is
what sets the streaming build's peak disk. At K=12 aa on GTEx that pool is ~1.93T pairs, so
two UInt32 per pair is ~15.4 TB on its own. The narrow format stores both as UInt16 and
halves it. Samples fit by construction (2,486 GTEx samples). Counts >= 0xffff go to an
exception list (UInt64 index + UInt32 value, 12 bytes each). Measured on the K=10 GTEx table
that is ~9e-6 of pairs (~9.1M, about 110 MB at GTEx scale) and 1.9e-7 on the K=12 small
cohort, so the escape is noise against 7.7 TB saved.

Narrow files start with the magic Int64 -2. Legacy files start with n_kmers >= 0, so the
reader takes either and a work dir from before this change still resumes. Only the disk
format changes, the in-memory `_PartialShard` stays UInt32. A shard whose max sample slot
exceeds 0xffff (a >65,535-sample cohort) falls back to the legacy layout. =#
const _PARTIAL_NARROW_MAGIC = Int64(-2)
const _PARTIAL_C16_ESC = UInt16(0xffff)

function _write_partial_shard(path::String, ps::_PartialShard)
    narrow = isempty(ps.pool_s) || maximum(ps.pool_s) <= typemax(UInt16)
    open(path, "w") do io
        narrow && write(io, _PARTIAL_NARROW_MAGIC)
        write(io, Int64(length(ps.kmers)))
        write(io, ps.kmers)
        write(io, ps.row_id)
        write(io, Int64(length(ps.offsets) - 1))  # n_rows
        write(io, ps.offsets)
        write(io, Int64(length(ps.pool_s)))  # n_pairs
        if !narrow
            write(io, ps.pool_s)
            write(io, ps.pool_c)
            return
        end
        np = length(ps.pool_s)
        buf = Vector{UInt16}(undef, np)
        @inbounds for i in 1:np; buf[i] = ps.pool_s[i] % UInt16; end
        write(io, buf)
        exc_i = UInt64[]; exc_v = UInt32[]
        @inbounds for i in 1:np
            c = ps.pool_c[i]
            if c >= _PARTIAL_C16_ESC
                buf[i] = _PARTIAL_C16_ESC
                push!(exc_i, UInt64(i)); push!(exc_v, c)
            else
                buf[i] = c % UInt16
            end
        end
        write(io, buf)
        write(io, Int64(length(exc_i)))
        write(io, exc_i)
        write(io, exc_v)
    end
end

function _read_partial_shard(path::String)::_PartialShard
    open(path, "r") do io
        first = read(io, Int64)
        narrow = first == _PARTIAL_NARROW_MAGIC
        nk = narrow ? read(io, Int64) : first
        kmers = Vector{UInt64}(undef, nk); read!(io, kmers)
        row_id = Vector{UInt32}(undef, nk); read!(io, row_id)
        nr = read(io, Int64)
        offsets = Vector{UInt64}(undef, nr + 1); read!(io, offsets)
        np = read(io, Int64)
        if !narrow
            pool_s = Vector{UInt32}(undef, np); read!(io, pool_s)
            pool_c = Vector{UInt32}(undef, np); read!(io, pool_c)
            return _PartialShard(kmers, row_id, offsets, pool_s, pool_c)
        end
        buf = Vector{UInt16}(undef, np)
        read!(io, buf); pool_s = UInt32.(buf)
        read!(io, buf); pool_c = UInt32.(buf)
        buf = nothing
        ne = read(io, Int64)
        exc_i = Vector{UInt64}(undef, ne); read!(io, exc_i)
        exc_v = Vector{UInt32}(undef, ne); read!(io, exc_v)
        @inbounds for j in 1:ne
            pool_c[exc_i[j]] == _PARTIAL_C16_ESC || error("$path: corrupt count exception at pair $(exc_i[j])")
            pool_c[exc_i[j]] = exc_v[j]
        end
        return _PartialShard(kmers, row_id, offsets, pool_s, pool_c)
    end
end

## Row pool builder (shared by wave build, merge and assembly) ##

# Accumulates distinct sparse rows. `intern!(pairs)` returns the row id for a
# slot-ascending `Vector{Tuple{UInt32,UInt32}}`, adding it to the pool on first sight.
mutable struct _RowPool
    offsets::Vector{UInt64}
    pool_s::Vector{UInt32}
    pool_c::Vector{UInt32}
    dedup::Dict{Vector{UInt32}, UInt32}
end
_RowPool() = _RowPool(UInt64[1], UInt32[], UInt32[], Dict{Vector{UInt32}, UInt32}())

function _intern_row!(rp::_RowPool, pairs::Vector{Tuple{UInt32, UInt32}})::UInt32
    key = Vector{UInt32}(undef, 2 * length(pairs))
    @inbounds for (t, (s, c)) in enumerate(pairs)
        key[2t - 1] = s
        key[2t] = c
    end
    return get!(rp.dedup, key) do
        @inbounds for (s, c) in pairs
            push!(rp.pool_s, s); push!(rp.pool_c, c)
        end
        push!(rp.offsets, UInt64(length(rp.pool_s) + 1))
        UInt32(length(rp.offsets) - 1)
    end
end

## Wave: count -> shard record files -> partial shards ##

# A transient counter failure (busy-NFS short read, samtools SIGPIPE, FIFO desync) used to
# kill a multi-day build outright; retry a few times before giving up and naming the sample.
function _count_with_retry(counter, path::AbstractString, slot::Integer, attempts::Int)::Dict{UInt64, UInt32}
    local err
    for a in 1:attempts
        try
            return counter(path)::Dict{UInt64, UInt32}
        catch e
            err = e
            a < attempts || break
            _log("  ! sample $slot count failed (attempt $a/$attempts): ", sprint(showerror, e))
            _log("  ! retrying $(basename(path)) in $(5a)s")
            sleep(5a)
        end
    end
    error("sample $slot ($path) failed to count after $attempts attempts; last error: " *
          sprint(showerror, err))
end

# Runs one wave: counts `paths` (global sample indices first_global .. first_global+len-1),
# streams records into per-shard files under `rec_dir`, then folds each shard into
# `out_dir/shard{p}.bin`. `counter(path)::Dict{UInt64,UInt32}`.
function _run_wave(paths::Vector{String}, first_global::Int, out_dir::String, rec_dir::String,
                   keep_shift::Int, P::Int, counter; count_attempts::Int = 4)
    mkpath(rec_dir)
    mask = UInt64(P - 1)
    flush_at = 1_000_000
    ios = IOStream[open(joinpath(rec_dir, "rec$(p).bin"), "w") for p in 0:(P - 1)]
    bufs = [_StreamRec[] for _ in 0:(P - 1)]
    try
        for (li, path) in enumerate(paths)
            slot = UInt32(first_global - 1 + li)
            t0 = time()
            ht = _count_with_retry(counter, path, slot, count_attempts)
            for (k, c) in ht
                p = _shard_of(k, keep_shift, mask) + 1
                push!(bufs[p], _StreamRec(k, slot, c))
                if length(bufs[p]) >= flush_at
                    write(ios[p], bufs[p]); empty!(bufs[p])
                end
            end
            _log("  sample $li/$(length(paths)) (global $slot): $(length(ht)) k-mers in ",
                 round(time() - t0; digits = 1), "s  ", basename(path))
            ht = nothing
        end
        for p in 1:P
            !isempty(bufs[p]) && write(ios[p], bufs[p])
        end
    finally
        for io in ios; close(io); end
    end

    rec_bytes = _dir_bytes(rec_dir)
    _log("  folding $P shard record files into partials ($(Base.format_bytes(rec_bytes)) of records)")
    for p in 0:(P - 1)
        rp = joinpath(rec_dir, "rec$(p).bin")
        _build_wave_shard(rp, joinpath(out_dir, "shard$(p).bin"))
        rm(rp; force = true)
    end
    open(joinpath(out_dir, "manifest"), "w") do io
        write(io, Int64(first_global)); write(io, Int64(length(paths)))
    end
    return rec_bytes
end

# Sort one shard's records by k-mer, collapse equal-k-mer runs into sparse rows, dedup
# rows, and write the shard's partial.
function _build_wave_shard(rec_path::String, out_path::String)
    n = Int(filesize(rec_path) ÷ sizeof(_StreamRec))
    recs = Vector{_StreamRec}(undef, n)
    if n > 0
        open(rec_path, "r") do io; read!(io, recs); end
        psort!(recs)
    end

    kmers = UInt64[]
    row_id = UInt32[]
    rp = _RowPool()

    i = 1
    while i <= n
        k = recs[i].kmer
        j = i
        pairs = Tuple{UInt32, UInt32}[]
        @inbounds while j <= n && recs[j].kmer == k
            push!(pairs, (recs[j].slot, recs[j].count))
            j += 1
        end
        sort!(pairs; by = first)
        push!(kmers, k)
        push!(row_id, _intern_row!(rp, pairs))
        i = j
    end
    _write_partial_shard(out_path, _PartialShard(kmers, row_id, rp.offsets, rp.pool_s, rp.pool_c))
end

## k-way merge of same-index shards from several partials ##

# Merge `shards` (all the same prefix shard, k-mers sorted within each) into one partial
# shard: for every k-mer, gather its (sample, count) pairs from every input holding it,
# order by sample, and dedup the combined row.
function _merge_partial_shards(shards::Vector{_PartialShard})::_PartialShard
    ns = length(shards)
    cur = ones(Int, ns)  # next k-mer index per shard
    lens = Int[length(s.kmers) for s in shards]

    out_kmers = UInt64[]
    out_row_id = UInt32[]
    rp = _RowPool()

    while true
        mink = typemax(UInt64); live = false
        @inbounds for si in 1:ns
            if cur[si] <= lens[si]
                live = true
                k = shards[si].kmers[cur[si]]
                k < mink && (mink = k)
            end
        end
        live || break

        pairs = Tuple{UInt32, UInt32}[]
        @inbounds for si in 1:ns
            (cur[si] <= lens[si] && shards[si].kmers[cur[si]] == mink) || continue
            sh = shards[si]
            r = sh.row_id[cur[si]]
            for t in sh.offsets[r]:(sh.offsets[r + 1] - 1)
                push!(pairs, (sh.pool_s[t], sh.pool_c[t]))
            end
            cur[si] += 1
        end
        sort!(pairs; by = first)
        push!(out_kmers, mink)
        push!(out_row_id, _intern_row!(rp, pairs))
    end
    return _PartialShard(out_kmers, out_row_id, rp.offsets, rp.pool_s, rp.pool_c)
end

# Read shard `p` from every partial dir in `dirs` and merge them. Shards are independent, so
# a merge node runs one of these per shard across all threads (see _run_merge_node!).
_merge_one_shard(dirs::Vector{String}, p::Int)::_PartialShard =
    _merge_partial_shards(_PartialShard[_read_partial_shard(joinpath(d, "shard$(p).bin")) for d in dirs])

## Hierarchical merge plan ##

# Bottom-up list of (output_name, input_names) super-partials, plus the final root list
# (<= fanin partials the assembly consumes directly).
function _merge_plan(names::Vector{String}, fanin::Int)
    plan = Tuple{String, Vector{String}}[]
    cur = names
    level = 0
    while length(cur) > fanin
        level += 1
        nxt = String[]
        for (g, grp) in enumerate(Iterators.partition(cur, fanin))
            out = "super_L$(level)_$(g)"
            push!(plan, (out, collect(String, grp)))
            push!(nxt, out)
        end
        cur = nxt
    end
    return plan, cur
end

## Adaptive prefix-tree reshard (fixes assembly's OOM on skewed shard sizes) ##

# Splits one on-disk partial shard into 2^extra_bits children by the next `extra_bits` bits
# below `remaining_bits` (bits not yet consumed by shard_bits or any prior split level).
# Reads the whole file once (same O(filesize) footprint as `_build_wave_shard`); callers may
# run this concurrently across the handful of ROOTS holding the same prefix (peak memory
# then bounded by root-count x one file's size, a few hundred GB at worst) but must never run
# it concurrently for two DIFFERENT prefixes -- the combined-across-roots size that triggers
# a split can be hundreds of GB, and stacking several different oversized shards at once is
# exactly the concurrency spike `_assemble_v4`'s own window caused in production (see
# `_resolve_leaf`, which parallelizes only across roots, never across prefixes). A k-mer's
# child index is monotonic in its sort position (the parent's `kmers` are
# already globally sorted and every k-mer in it shares identical already-consumed bits), so
# a single linear scan finds each child's contiguous index range in order, and each child is
# built and written before moving to the next, bounding peak extra memory to the parent plus
# its single largest child rather than a multiple of the parent.
#
# Rows are copied flat, one row per k-mer, deliberately WITHOUT `_RowPool`/`_intern_row!`
# cross-k-mer dedup -- same no-dedup tradeoff `_merge_shard_rows` already makes at assembly,
# for the same measured reason (see ARCHITECTURE.md's CountsLayer section: rows longer than
# 4 pairs dedup at ~1.0x in practice, essentially unique, so the dedup dictionary costs about
# as many bytes as it saves). It also mattered in practice: an earlier Dict-based version of
# this function took ~10.5 minutes on a real 43.8GB production shard (552M rows); this flat-
# copy version does the equivalent work as a handful of `copyto!` calls per child.
function _split_shard_file!(src_path::String, dest_dir::String, out_prefix::String,
                            remaining_bits::Int, extra_bits::Int)
    ps = _read_partial_shard(src_path)
    nk = length(ps.kmers)
    child_shift = remaining_bits - extra_bits
    nchildren = 1 << extra_bits
    child_mask = UInt64(nchildren - 1)
    cidx(i) = Int((ps.kmers[i] >> child_shift) & child_mask)

    lo = 1
    for c in 0:(nchildren - 1)
        hi = lo
        while hi <= nk && cidx(hi) == c
            hi += 1
        end
        n = hi - lo
        kmers = Vector{UInt64}(undef, n)
        row_id = Vector{UInt32}(undef, n)
        offsets = Vector{UInt64}(undef, n + 1)
        offsets[1] = 1
        @inbounds for (j, i) in enumerate(lo:(hi - 1))
            kmers[j] = ps.kmers[i]
            row_id[j] = UInt32(j)
            r = ps.row_id[i]
            offsets[j + 1] = offsets[j] + (ps.offsets[r + 1] - ps.offsets[r])
        end
        pool_s = Vector{UInt32}(undef, offsets[end] - 1)
        pool_c = Vector{UInt32}(undef, offsets[end] - 1)
        @inbounds for (j, i) in enumerate(lo:(hi - 1))
            r = ps.row_id[i]
            a = ps.offsets[r]; len = ps.offsets[r + 1] - a
            dst = offsets[j]
            copyto!(pool_s, dst, ps.pool_s, a, len)
            copyto!(pool_c, dst, ps.pool_c, a, len)
        end
        _write_partial_shard(joinpath(dest_dir, "$(out_prefix)_$(c).bin"),
                             _PartialShard(kmers, row_id, offsets, pool_s, pool_c))
        lo = hi
    end
    return nothing
end

# Decide/split/recurse for one leaf candidate (a base shard or an already-split child).
# `paths` is that prefix's already-resolved file per root. Splits only when the COMBINED
# size across every root exceeds `leaf_bytes` (never per-root -- assembly reads matching
# leaf indices from every root at once, so an asymmetric split would misalign them) and
# enough bits remain. `work/reshard/<rootbasename>/` is shared by the whole tree for that
# root, so a retry only overwrites this prefix's own child files, never wipes the directory
# (that would destroy already-completed sibling splits).
#
# The per-root splits run concurrently (`Threads.@spawn`, one task per root): each root's
# file is independent (separate source file, separate `dest` directory), so there is no
# shared mutable state to race on, and peak memory stays bounded by root-count x one file's
# size (a handful of roots, a few hundred GB worst case) rather than growing with however
# many *different* prefixes happen to be oversized -- this call is only ever reached for one
# prefix at a time (`_reshard_leaves`'s outer loop and this function's own recursion into
# children are both strictly sequential), which is the property that actually matters for
# memory safety.
function _resolve_leaf(root_dirs::Vector{String}, paths::Vector{String}, prefix::String,
                       remaining_bits::Int, work::String, leaf_bytes::Int, extra_bits::Int,
                       _done, _mark)::Vector{Vector{String}}
    total = sum(filesize, paths)
    (total <= leaf_bytes || remaining_bits < extra_bits) && return [paths]

    marker = "reshard_$(prefix)"
    if !_done(marker)
        _log("  reshard: splitting $(prefix) ($(Base.format_bytes(total)) combined) into ",
             "$(1 << extra_bits) children across $(length(root_dirs)) roots")
        t0 = time()
        tasks = [Threads.@spawn begin
                    dest = joinpath(work, "reshard", basename(d))
                    mkpath(dest)
                    _split_shard_file!(src, dest, prefix, remaining_bits, extra_bits)
                 end for (d, src) in zip(root_dirs, paths)]
        foreach(fetch, tasks)
        _mark(marker)
        _log("    $(prefix) split in ", round(time() - t0; digits = 1), "s")
    end

    out = Vector{Vector{String}}()
    for c in 0:((1 << extra_bits) - 1)
        cprefix = "$(prefix)_$(c)"
        cpaths = [joinpath(work, "reshard", basename(d), "$(cprefix).bin") for d in root_dirs]
        append!(out, _resolve_leaf(root_dirs, cpaths, cprefix, remaining_bits - extra_bits,
                                   work, leaf_bytes, extra_bits, _done, _mark))
    end
    return out
end

# Top-level entry: walks every base shard 0:(P-1) in order, resolving each to one or more
# leaves. Ascending order (both here and in `_resolve_leaf`'s child loop) preserves the
# global k-mer ordering `_assemble_v4` depends on, since a parent's k-mers are already
# sorted and a split only ever refines, never reorders, that ordering.
function _reshard_leaves(root_dirs::Vector{String}, P::Int, keep_shift::Int, work::String,
                         leaf_bytes::Int, extra_bits::Int, _done, _mark)::Vector{Vector{String}}
    _log("reshard: checking $P base shards (leaf target $(Base.format_bytes(leaf_bytes)), " *
         "$(1 << extra_bits)-way fanout)")
    leaves = Vector{Vector{String}}()
    for p in 0:(P - 1)
        prefix = "shard$(p)"
        paths = [joinpath(d, "$(prefix).bin") for d in root_dirs]
        append!(leaves, _resolve_leaf(root_dirs, paths, prefix, keep_shift, work,
                                      leaf_bytes, extra_bits, _done, _mark))
        (p + 1) % 32 == 0 && _log("  reshard: $(p + 1)/$P base shards checked, " *
                                  "$(length(leaves)) leaves so far")
    end
    _log("reshard: done, $(length(leaves)) leaves from $P base shards")
    return leaves
end

## Assembly of the final .kct ##

const _V4_BLOCK_KMERS = 64  # k-mers per packed block in the CountsLayer

# Merge one leaf's rows from every root partial into per-k-mer rows, WITHOUT deduplication,
# in flat CSR form (row i is `flat_s[off[i]:off[i+1]-1]` / `flat_c[...]`, no Vector-per-k-mer).
# Roots cover disjoint, ascending sample ranges and are already sample-sorted, so a k-mer's
# row is just its per-root pieces concatenated (see ARCHITECTURE.md for why no dedup here,
# unlike `_merge_partial_shards` at the wave/super level). `paths` is one already-resolved
# file per root for this leaf (a base `shard{p}.bin`, or a reshard child under `work/reshard`
# if that shard was oversized -- see `_reshard_leaves`).
function _merge_shard_rows(paths::Vector{String})
    shards = _PartialShard[_read_partial_shard(p) for p in paths]
    ns = length(shards)
    cur = ones(Int, ns)
    lens = Int[length(s.kmers) for s in shards]

    out_k = UInt64[]
    off = Int[1]
    flat_s = UInt32[]
    flat_c = UInt32[]

    while true
        mink = typemax(UInt64); live = false
        @inbounds for si in 1:ns
            if cur[si] <= lens[si]
                live = true
                k = shards[si].kmers[cur[si]]
                k < mink && (mink = k)
            end
        end
        live || break

        @inbounds for si in 1:ns
            (cur[si] <= lens[si] && shards[si].kmers[cur[si]] == mink) || continue
            sh = shards[si]
            r = sh.row_id[cur[si]]
            rng = sh.offsets[r]:(sh.offsets[r + 1] - 1)
            append!(flat_s, @view sh.pool_s[rng])
            append!(flat_c, @view sh.pool_c[rng])
            cur[si] += 1
        end
        push!(out_k, mink)
        push!(off, length(flat_s) + 1)
    end
    return out_k, off, flat_s, flat_c
end

#= Incremental DeltaArray encoder for assembly. The old path appended every k-mer into one
Vector{UInt64} and only then delta-encoded it: 8 bytes per k-mer, up to ~16 while append!
grows the vector, plus the 4-byte deltas built next to it. That is why the K = 10 GTEx assembly
(62.5 B k-mers) needed an 1,800 GB job, and a K = 12 cut-and-pad table with more k-mers would
not fit on oni at all. Encoding as k-mers stream out of the leaves keeps only the deltas
(4 bytes each, ~8 while growing) plus the checkpoints. Output is identical to NArrays'
_delta_encode on the same sorted sequence (checked in the test suite). =#
mutable struct _DeltaBuilder
    checkpoints::Vector{UInt64}
    deltas::Vector{UInt32}
    regular_cp_idx::Vector{UInt64}
    inter::Int
    prev::UInt64
    n::Int
end
_DeltaBuilder(inter::Int) = _DeltaBuilder(UInt64[], UInt32[], UInt64[], inter, UInt64(0), 0)

function Base.push!(b::_DeltaBuilder, v::UInt64)
    b.n += 1
    b.n > 1 && v <= b.prev && error("k-mers out of order at $(b.n): $v after $(b.prev)")
    is_regular = (b.n - 1) % b.inter == 0
    if is_regular || (b.n > 1 && v - b.prev > typemax(UInt32))
        is_regular && push!(b.regular_cp_idx, length(b.checkpoints) + 1)
        push!(b.checkpoints, v)
        push!(b.deltas, UInt32(0))
    else
        push!(b.deltas, UInt32(v - b.prev))
    end
    b.prev = v
    return b
end

Base.length(b::_DeltaBuilder) = b.n
_finish(b::_DeltaBuilder) = DeltaArray(b.checkpoints, b.deltas, b.regular_cp_idx, b.inter)

# Transcode the merged root partials into one V4.0 table, packing each k-mer's row inline.
# The blob streams to a temp file next to `out_path` rather than staying resident (it can be
# 1-2 TB at full cohort, see ARCHITECTURE.md); RAM holds the k-mer array, block index and one
# merge window. Per-leaf merges run in a bounded parallel window; the pack loop is serial.
# `leaves` is the ordered, already-resolved leaf list from `_reshard_leaves` (ascending
# global k-mer order): each entry is one already-resolved file path per root. A leaf is
# either an original `shard{p}.bin` left alone, or a reshard child under `work/reshard`
# for a shard that was oversized and split further -- see ARCHITECTURE.md.
function _assemble_v4(leaves::Vector{Vector{String}}, kk::Int, Ab::Type, idx_prefix::Int,
                     total_samples::Int, out_path::String; nk_hint::Int = 0)
    B = _V4_BLOCK_KMERS
    L = length(leaves)
    blobtmp = string(out_path, ".blob.tmp")
    rm(blobtmp; force = true)

    kbuild = _DeltaBuilder(DEFAULT_CHECKPOINT_INTERVAL)
    #= Growing `deltas` by push! reallocates and copies, so near the end the old and new
    buffers coexist: at K=12 aa GTEx that's ~512 GB live plus a ~768 GB copy target, right at
    the job's memory limit. `nk_hint` is the k-mer count summed over roots, an overcount
    (a k-mer in several roots counts several times), so this reserves more than it fills.
    Linux only commits pages that get written, so the slack costs address space, not RSS
    (measured, see PAPER_TODO D6). =#
    nk_hint > 0 && sizehint!(kbuild.deltas, nk_hint)
    block_ptr = UInt64[]
    npairs = 0
    blob_len = 0
    W = clamp(Threads.nthreads(), 1, 8)

    open(blobtmp, "w") do bio
        blobbuf = UInt8[]
        # current partial block, flat: row j is pend_s[pend_off[j]:pend_off[j+1]-1]
        pend_s = UInt32[]; pend_c = UInt32[]; pend_off = Int[1]
        function flush_block!()
            n = length(pend_off) - 1
            n == 0 && return
            push!(block_ptr, UInt64(blob_len))
            before = length(blobbuf)
            _pack_block!(blobbuf, pend_s, pend_c, pend_off, n)
            blob_len += length(blobbuf) - before
            empty!(pend_s); empty!(pend_c); resize!(pend_off, 1)  # pend_off[1] stays 1
            if length(blobbuf) >= (1 << 26)
                write(bio, blobbuf); empty!(blobbuf)
            end
        end

        nleaves_done = 0
        for chunk in Iterators.partition(leaves, W)
            tasks = [Threads.@spawn _merge_shard_rows(leafpaths) for leafpaths in chunk]
            for t in tasks
                ks, roff, fs, fc = fetch(t)
                foreach(k -> push!(kbuild, k), ks)
                @inbounds for i in eachindex(ks)
                    a = roff[i]; b = roff[i + 1] - 1
                    npairs += (b - a + 1)
                    append!(pend_s, @view fs[a:b])
                    append!(pend_c, @view fc[a:b])
                    push!(pend_off, length(pend_s) + 1)
                    length(pend_off) - 1 == B && flush_block!()
                end
            end
            nleaves_done += length(chunk)
            _log("  assemble: $nleaves_done/$L leaves packed, $(length(kbuild)) k-mers, " *
                 "$(Base.format_bytes(blob_len)) blob, RSS $(round(_rss_gb(); digits = 1)) GB")
        end
        flush_block!()
        isempty(blobbuf) || write(bio, blobbuf)
    end

    n_kmers = length(kbuild)
    seqs = _finish(kbuild)  # encoded as the leaves streamed out, globally sorted
    kbuild = nothing
    kl = KmerLayer{kk, Ab, UInt64, UInt32}(seqs, _empty_kmer_idx(kk, Ab(), UInt64; prefix_size = idx_prefix))
    cl = CountsLayer(B, block_ptr, UInt8[], n_kmers, total_samples)  # blob stays on disk
    kct = KCT(kl, cl)
    compute_index!(kct.kmer; prefix_size = idx_prefix)
    _log("  assembled: $n_kmers k-mers, $npairs pairs, blob $(Base.format_bytes(blob_len)), " *
         "block_ptr $(Base.format_bytes(8 * length(block_ptr)))")

    open(out_path, "w") do io
        write(io, 4.0)
        _write_header_and_kmers(io, kct)
        _write_counts_header(io, total_samples, n_kmers, B, length(block_ptr), blob_len, block_ptr)
        _move_blob_tail!(io, blobtmp, blob_len)
        _write_biotype(io, nothing)
    end
    rm(blobtmp; force = true)
    return out_path
end

#= The blob has to land after the k-mer section and block index, and neither size is known
until the last leaf is packed, so it goes to `blobtmp` first. A forward copy then holds the
blob twice (~4 TB extra at K=12 aa GTEx). This copies it back to front instead, writing each
64 MB chunk at its final offset in `io` (the file is sparse until filled) and truncating
`blobtmp` behind it, so the extra disk is one chunk. A crash mid-move loses the blob, which
is fine: assembly isn't resumable anyway and reruns from the leaves, which `work` still
holds. Leaves `io` positioned right after the blob. =#
function _move_blob_tail!(io::IO, blobtmp::String, blob_len::Integer; chunk::Int = 1 << 26)
    pos0 = position(io)
    open(blobtmp, "r+") do bin
        filesize(blobtmp) == blob_len || error("$blobtmp is $(filesize(blobtmp)) bytes, expected $blob_len")
        buf = Vector{UInt8}(undef, chunk)
        hi = Int(blob_len)
        nmoved = 0
        while hi > 0
            lo = max(0, hi - chunk)
            n = hi - lo
            seek(bin, lo)
            read!(bin, view(buf, 1:n))
            seek(io, pos0 + lo)
            write(io, view(buf, 1:n))
            flush(io)  # must be in the output before the source bytes go away
            truncate(bin, lo)
            hi = lo
            nmoved += 1
            nmoved % 4096 == 0 && _log("  blob move: $(Base.format_bytes(blob_len - hi)) of " *
                                       "$(Base.format_bytes(blob_len)) moved, RSS $(round(_rss_gb(); digits = 1)) GB")
        end
    end
    seek(io, pos0 + blob_len)
    return nothing
end

## Seed adapter: an existing .kct -> a level-0 partial ##

# Sparse (global_slot, count) pairs for k-mer i: decode its inline row. Decodes the whole
# block, so a full seed walk is O(block_size) per k-mer; acceptable for the occasional seed fold.
function _seed_row_pairs(cl::CountsLayer, i::Int, base::Int)
    s, c = _decode_row(cl, i)
    pairs = Tuple{UInt32, UInt32}[]
    @inbounds for j in eachindex(s)
        push!(pairs, (UInt32(base + s[j]), c[j]))
    end
    return pairs
end

# Turn an existing .kct into a partial table under `out_dir`, its samples numbered from
# `first_global_sample`. Returns the number of samples it contributes.
function stream_kct_as_shards(kct_path::String, first_global_sample::Int, out_dir::String,
                              keep_shift::Int, P::Int)
    kct = load_kct(kct_path)
    counts = kct.counts
    counts === nothing && error("seed KCT has no counts layer: $kct_path")
    n_src = Int(counts.n_samples.x)
    base = first_global_sample - 1

    kmers = [UInt64[] for _ in 0:(P - 1)]
    rids = [UInt32[] for _ in 0:(P - 1)]
    rps = [_RowPool() for _ in 0:(P - 1)]
    mask = UInt64(P - 1)

    it = iterate(kct.kmer.seqs); i = 1
    while it !== nothing
        kbits, st = it
        p = _shard_of(kbits, keep_shift, mask) + 1
        pairs = _seed_row_pairs(counts, i, base)
        sort!(pairs; by = first)
        push!(kmers[p], kbits)
        push!(rids[p], _intern_row!(rps[p], pairs))
        it = iterate(kct.kmer.seqs, st); i += 1
    end

    for p in 1:P
        _write_partial_shard(joinpath(out_dir, "shard$(p - 1).bin"),
                             _PartialShard(kmers[p], rids[p], rps[p].offsets, rps[p].pool_s, rps[p].pool_c))
    end
    open(joinpath(out_dir, "manifest"), "w") do io
        write(io, Int64(first_global_sample)); write(io, Int64(n_src))
    end
    return n_src
end

## Single-sample patch (post-assembly correction, no resharding, no new sample slot) ##

# Replace one sample's contribution in an already-assembled .kct with a fresh, correct
# recount, writing a NEW file at `out_path` (never touches `kct_path`). Built for the
# one-off GTEx DS rebuild anomaly (slot 289/SRR614719 undercounted by a transient build-time
# glitch, 28.1M vs the correct 335.5M k-mers -- see PAPER_TODO.md C11): a recount that was
# never truncated, like the one `cram_hash` produces fresh, is a strict superset of whatever
# a truncated build-time count wrote for the same sample, so this never needs to "merge"
# counts for a single (kmer, sample) pair, only replace them outright.
#
# Streams the existing table once in its existing sorted order against the sorted corrected
# recount -- a union merge, the same shape as `_merge_shard_rows`'s k-way merge -- and
# repacks every row exactly as `_assemble_v4` does. A k-mer only the OLD (wrong) count
# supported and the new recount does not is dropped entirely if that empties its row (should
# not happen in practice given the superset property above, but costs nothing to handle).
function patch_sample_counts(kct_path::String, sample_slot::Integer,
                             corrected_counts::Dict{UInt64, UInt32}, out_path::String)
    kct = load_kct(kct_path)
    kct.counts === nothing && error("source KCT has no counts layer: $kct_path")
    kct.biotype === nothing ||
        error("source KCT has a biotype layer: patching reorders/add k-mers and would " *
              "desync a per-k-mer-index biotype annotation -- not supported here")
    return _patch_sample_counts(kct, Int(sample_slot), corrected_counts, out_path)
end

function _patch_sample_counts(kct::KCT{K, Ab, CountsLayer, Nothing, C, D}, sample_slot::Int,
                              corrected_counts::Dict{UInt64, UInt32},
                              out_path::String) where {K, Ab <: Alphabet, C <: Unsigned, D <: Unsigned}
    cl = kct.counts
    n_samples = Int(cl.n_samples.x)
    1 <= sample_slot <= n_samples || error("sample_slot=$sample_slot out of range 1:$n_samples")
    slot = UInt32(sample_slot)
    idx_prefix = idx_prefix_size(kct)
    n_existing = length(kct)

    corrected_kmers = sort!(collect(keys(corrected_counts)))
    nc = length(corrected_kmers)
    _log("patch_sample_counts: sample $sample_slot, $n_existing existing k-mers, " *
         "$nc corrected k-mers")

    B = _V4_BLOCK_KMERS
    blobtmp = string(out_path, ".blob.tmp")
    rm(blobtmp; force = true)

    all_kmers = UInt64[]
    block_ptr = UInt64[]
    npairs = 0
    blob_len = 0
    n_scrubbed = 0
    n_touched = 0
    n_new = 0

    open(blobtmp, "w") do bio
        blobbuf = UInt8[]
        pend_s = UInt32[]; pend_c = UInt32[]; pend_off = Int[1]
        function flush_block!()
            n = length(pend_off) - 1
            n == 0 && return
            push!(block_ptr, UInt64(blob_len))
            before = length(blobbuf)
            _pack_block!(blobbuf, pend_s, pend_c, pend_off, n)
            blob_len += length(blobbuf) - before
            empty!(pend_s); empty!(pend_c); resize!(pend_off, 1)
            if length(blobbuf) >= (1 << 26)
                write(bio, blobbuf); empty!(blobbuf)
            end
        end
        function emit_row!(kmer::UInt64, s::Vector{UInt32}, c::Vector{UInt32})
            push!(all_kmers, kmer)
            npairs += length(s)
            append!(pend_s, s); append!(pend_c, c)
            push!(pend_off, length(pend_s) + 1)
            length(pend_off) - 1 == B && flush_block!()
        end

        # Driven by `_each_row`, which decodes each 64-k-mer block exactly once (unlike
        # `_decode_row`, which is fine for scattered point lookups but would redecode the
        # same block up to 64x over in a full-table walk like this one -- the table's own
        # doc comment on `_decode_row` calls this out explicitly). `kseqs`/`st` track the
        # matching k-mer VALUE in lockstep, since `_each_row` only hands back `(i, s, c)`;
        # both advance by exactly one every callback, in the same 1..n_existing order.
        ci = 1
        kseqs = kct.kmer.seqs
        it = iterate(kseqs)
        _each_row(cl) do i, s, c
            kbits, st = it
            # corrected-only k-mers sorting strictly before this one: brand-new rows.
            while ci <= nc && corrected_kmers[ci] < kbits
                n_new += 1
                kc = corrected_kmers[ci]
                emit_row!(kc, UInt32[slot], UInt32[corrected_counts[kc]])
                ci += 1
            end
            if ci <= nc && corrected_kmers[ci] == kbits
                n_touched += 1
                j = findfirst(==(slot), s)
                s2 = j === nothing ? copy(s) : deleteat!(copy(s), j)
                c2 = j === nothing ? copy(c) : deleteat!(copy(c), j)
                ins = searchsortedfirst(s2, slot)
                insert!(s2, ins, slot); insert!(c2, ins, corrected_counts[kbits])
                emit_row!(kbits, s2, c2)
                ci += 1
            else
                j = findfirst(==(slot), s)
                if j === nothing
                    emit_row!(kbits, s, c)
                else
                    n_scrubbed += 1
                    s2 = deleteat!(copy(s), j); c2 = deleteat!(copy(c), j)
                    isempty(s2) || emit_row!(kbits, s2, c2)
                end
            end
            it = iterate(kseqs, st)

            if i % 5_000_000 == 0
                _log("  patch: $i/$n_existing existing k-mers walked, $(length(all_kmers)) " *
                     "emitted, $(Base.format_bytes(blob_len)) blob, RSS $(round(_rss_gb(); digits=1)) GB")
            end
        end
        # corrected-only k-mers sorting after every existing k-mer.
        while ci <= nc
            n_new += 1
            kc = corrected_kmers[ci]
            emit_row!(kc, UInt32[slot], UInt32[corrected_counts[kc]])
            ci += 1
        end
        flush_block!()
        isempty(blobbuf) || write(bio, blobbuf)
    end

    n_kmers = length(all_kmers)
    seqs = DeltaArray{C, D}(all_kmers, DEFAULT_CHECKPOINT_INTERVAL)  # already globally sorted
    all_kmers = UInt64[]
    kl = KmerLayer{K, Ab, C, D}(seqs, _empty_kmer_idx(K, Ab(), C; prefix_size = idx_prefix))
    cl_new = CountsLayer(B, block_ptr, UInt8[], n_kmers, n_samples)
    kct_new = KCT(kl, cl_new)
    compute_index!(kct_new.kmer; prefix_size = idx_prefix)
    _log("patch_sample_counts: done. $n_kmers k-mers ($(n_kmers - n_existing >= 0 ? "+" : "")" *
         "$(n_kmers - n_existing) net), $n_scrubbed rows had sample $sample_slot scrubbed, " *
         "$n_touched rows updated in place, $n_new brand-new rows, " *
         "blob $(Base.format_bytes(blob_len))")

    open(out_path, "w") do io
        write(io, 4.0)
        _write_header_and_kmers(io, kct_new)
        _write_counts_header(io, n_samples, n_kmers, B, length(block_ptr), blob_len, block_ptr)
        open(blobtmp, "r") do bin
            buf = Vector{UInt8}(undef, 1 << 26)
            while !eof(bin)
                n = readbytes!(bin, buf)
                write(io, view(buf, 1:n))
            end
        end
        _write_biotype(io, nothing)
    end
    rm(blobtmp; force = true)
    return out_path
end

## Free-space guard ##

#= The August GTEx DNA run wedged in an uninterruptible NFS write() when the export filled,
and a full disk mid-merge can leave a half-written shard behind. So before each step that
writes a lot, check the out_dir filesystem has room for that step plus `margin`, and if not,
log and sleep until someone frees space instead of failing. A paused job costs wall time, a
wedged one costs the run. Needs are estimates: a merge writes about what its inputs hold, and
assembly writes the blob (~2.1 bytes/pair measured on the K=10 GTEx table, 2.5 used here)
plus ~4 bytes per k-mer (deltas), with k-mers summed over roots so it overcounts. =#
_free_bytes(path::String) = Int(Base.Filesystem.diskstat(path).available)
_dir_bytes(d::String) = isdir(d) ? sum((filesize(joinpath(d, f)) for f in readdir(d)); init = 0) : 0

function _await_disk(path::String, need::Integer, what::String; poll_s::Int = 900)
    need <= 0 && return nothing
    waited = false
    while (f = _free_bytes(path)) < need
        _log("  DISK: $what needs $(Base.format_bytes(need)) free on $path, only " *
             "$(Base.format_bytes(f)) available. Waiting $(poll_s ÷ 60) min, free space to continue.")
        waited = true
        sleep(poll_s)
    end
    waited && _log("  DISK: $(Base.format_bytes(_free_bytes(path))) free now, continuing $what")
    return nothing
end

# (n_kmers, n_pairs) from a partial shard's headers, without reading its arrays.
function _partial_sizes(path::String)
    open(path, "r") do io
        nk = read(io, Int64)
        narrow = nk == _PARTIAL_NARROW_MAGIC
        narrow && (nk = read(io, Int64))
        skip(io, nk * (sizeof(UInt64) + sizeof(UInt32)))
        nr = read(io, Int64); skip(io, (nr + 1) * sizeof(UInt64))
        return nk, read(io, Int64)
    end
end

## Pool-size tripwire ##

function _partial_pool_bytes(dir::String, P::Int)
    total = 0
    for p in 0:(P - 1)
        open(joinpath(dir, "shard$(p).bin"), "r") do io
            nk = read(io, Int64)
            nk == _PARTIAL_NARROW_MAGIC && (nk = read(io, Int64))
            skip(io, nk * sizeof(UInt64) + nk * sizeof(UInt32))  # kmers + row_id
            nr = read(io, Int64); skip(io, (nr + 1) * sizeof(UInt64))  # offsets
            np = read(io, Int64)
            total += 2 * np * sizeof(UInt32)  # pool_s + pool_c
        end
    end
    return total
end

## Top-level ##

const _STREAM_TMP_DEFAULT = get(ENV, "NEOKCT_STREAM_TMP", "/scratch")

"""
    build_kct_streaming(sample_paths; out_dir, kwargs...) -> out_path

Build a V4.0 (`CountsLayer`) `.kct` from `sample_paths` by a k-way streaming merge.
Samples are counted in waves, bucketed by k-mer prefix into shard files, folded into
prefix-sharded partial tables, then combined by a hierarchical k-way merge and assembled
into one `.kct`. Nothing dense is ever resident, so peak memory stays flat across the
whole build. Each partial is merged upward as soon as `merge_fanin` siblings are ready
rather than after the last wave, so partial storage on disk stays at O(`merge_fanin`)
waves instead of growing with the cohort.

# Keyword Arguments
- `K`: nucleotide k-mer length. `chunks` is the read-chunk size passed to the counter
- `translate`: `false` (default) builds a DNA table, `KCT{K, DNAAlphabet{2}}`. `true` builds
  a translated table, `KCT{K÷3, AAAlphabet}`. The default `counter` picks up this flag
- `double_strand`: `false` (default, matches every table built before this option existed).
  `true` also counts each read's reverse complement, so a library where the sequenced strand
  isn't consistently the transcript's coding strand doesn't lose that coverage -- see
  `jello_superthreaded_hash`'s docstring for why this has to happen during counting rather
  than as a query-time trick once `translate=true` has folded reads down to AA k-mers. The
  default `counter` picks up this flag; a custom `counter` must honour it itself
- `pad`, `min_res`: cut-and-pad translated windows (needs `translate=true`), see
  `jello_superthreaded_hash`. Makes prefix queries exact read-occurrence counts and allows
  K up to 36 nt. The default `counter` picks these up
- `wave_size`: samples counted and folded per wave
- `shard_bits`: gives `2^shard_bits` prefix shards. `merge_fanin` partials are combined per
  hierarchical merge step
- `reshard_leaf_bytes`: after merges finish, any shard whose size summed across every root
  partial exceeds this is split further by `reshard_extra_bits` additional k-mer prefix
  bits (recursively, until under this target), so assembly's concurrency window never has
  to hold several oversized shards at once. Real amino-acid/nucleotide frequency is far
  from uniform, so a handful of prefixes can otherwise dominate -- see ARCHITECTURE.md
- `reshard_extra_bits`: fan-out per reshard split level (`2^reshard_extra_bits` children).
  Fixes the bit position every split level consumes, so it must stay the same across a
  resumed run against the same `out_dir` (checked and erred on mismatch, not silently
  honoured)
- `counter`: `counter(path)::Dict{UInt64,UInt32}` produces one sample's k-mer counts
  (default `jello_superthreaded_hash`, honouring `translate`)
- `seeds`: `Vector{Tuple{String,Int}}` of `(existing .kct, first_global_sample)` folded in as
  level-0 partials. Sample ranges must be disjoint and precede the waves
- `resume`: skips any wave or partial that already has a `DONE` marker
- `count_attempts`: retries per sample when `counter` throws (transient samtools/FIFO/NFS
  failures), with a linear backoff, before aborting and naming the sample
- `min_free_gb`: headroom kept on `out_dir`'s filesystem. Each wave, merge and the assembly
  wait (logging every 15 min) until free space covers their estimated output plus this
- `max_pool_gb`: accepted for CLI compatibility, no longer used (the build logs a size
  projection after the first super-partial but never aborts)
"""
function build_kct_streaming(sample_paths::Vector{String};
        K::Int = 30, chunks::Int = 500_000, translate::Bool = false, double_strand::Bool = false,
        pad::Bool = false, min_res::Int = 8,
        wave_size::Int = 100, shard_bits::Int = 8,
        merge_fanin::Int = 5, max_pool_gb::Real = 400, idx_prefix::Int = -1,
        reshard_leaf_bytes::Int = 4_000_000_000, reshard_extra_bits::Int = 4,
        count_attempts::Int = 4, min_free_gb::Real = 50,
        tmp_dir::String = _STREAM_TMP_DEFAULT, out_dir::String,
        out_path::String = joinpath(out_dir, "neokct_v4.kct"),
        counter = (p -> jello_superthreaded_hash(p, K, chunks; translate = translate, double_strand = double_strand,
                                                 pad = pad, min_res = min_res)),
        seeds::Vector{Tuple{String, Int}} = Tuple{String, Int}[],
        resume::Bool = true)

    Ab = translate ? AAAlphabet : DNAAlphabet{2}
    bps = translate ? _STREAM_BITS_PER_AA : _STREAM_BITS_PER_NT
    kk = translate ? K ÷ 3 : K
    # `idx_prefix` < 0 means derive it via the same alphabet-dispatched default `load_kct`
    # uses when reloading this table -- they must agree, or a fresh load silently allocates
    # the wrong bucket count (see `_default_idx_prefix_size` in KCTLayers.jl).
    idx_prefix < 0 && (idx_prefix = _default_idx_prefix_size(kk, Ab()))
    P = 1 << shard_bits
    keep_shift = bps * kk - shard_bits
    keep_shift >= 1 || error("shard_bits=$shard_bits too large for a $(kk)-symbol k-mer")
    merge_fanin >= 2 || error("merge_fanin must be >= 2")

    work = joinpath(out_dir, "streaming_work"); mkpath(work)
    margin = round(Int, min_free_gb * 1e9)
    rec_root = joinpath(tmp_dir, "neokct_stream_rec_$(getpid())"); mkpath(rec_root)

    # DONE markers live in a flat dir, not inside each partial's folder: merges now run
    # interleaved with the waves, so a wave's folder is deleted long before the build ends
    # and a marker kept inside it would be lost. Resume must still see that wave as done.
    mark_dir = joinpath(work, "_markers"); mkpath(mark_dir)
    _done(name) = isfile(joinpath(mark_dir, name))
    _mark(name) = touch(joinpath(mark_dir, name))
    covered = Dict{String, Int}()  # partial name -> samples it covers

    try
        # seeds
        seed_names = String[]
        n_prefix = 0
        for (si, (path, first_s)) in enumerate(seeds)
            name = "seed_$(si)"; d = joinpath(work, name)
            first_s == n_prefix + 1 || error("seed $si first_global_sample=$first_s, expected $(n_prefix + 1)")
            if !(resume && _done(name))
                rm(d; recursive = true, force = true); mkpath(d)
                nS = stream_kct_as_shards(path, first_s, d, keep_shift, P)
                write(joinpath(mark_dir, name * ".n"), string(nS))  # survives the folder
                _mark(name)
            end
            nfile = joinpath(mark_dir, name * ".n")
            isfile(nfile) || error("seed $si is marked done but $nfile is missing; delete $work and rebuild")
            nS = parse(Int, read(nfile, String))
            covered[name] = nS; n_prefix += nS
            push!(seed_names, name)
            _log("seed $si: $path -> $nS samples (global $(first_s)..$(n_prefix))")
        end

        # waves
        n_waves = cld(length(sample_paths), wave_size)
        wave_names = ["wave_$(w)" for w in 1:n_waves]
        final_total = n_prefix + length(sample_paths)

        # The whole merge tree over every leaf, known up front: node names, inputs and
        # sample coverage are a pure function of the seed + wave list, so `covered` can be
        # filled in now and never depends on a folder that a merge has since deleted.
        for w in 1:n_waves
            lo = (w - 1) * wave_size + 1
            hi = min(w * wave_size, length(sample_paths))
            covered["wave_$(w)"] = hi - lo + 1
        end
        plan, roots = _merge_plan(vcat(seed_names, wave_names), merge_fanin)
        for (out, ins) in plan
            covered[out] = sum(covered[i] for i in ins)
        end

        # Run one merge node: k-way merge its inputs per shard, mark it, delete the inputs.
        tripped = Ref(false)
        function _run_merge_node!(out::String, ins::Vector{String})
            od = joinpath(work, out)
            rm(od; recursive = true, force = true)
            _await_disk(work, sum(i -> _dir_bytes(joinpath(work, i)), ins) + margin, "merge $out")
            mkpath(od)
            _log("  merge $out <- [$(join(ins, ", "))]  ($(covered[out]) samples), $(Threads.nthreads()) threads")
            in_dirs = String[joinpath(work, i) for i in ins]
            ndone = Threads.Atomic{Int}(0)
            Threads.@threads for p in 0:(P - 1)
                _write_partial_shard(joinpath(od, "shard$(p).bin"), _merge_one_shard(in_dirs, p))
                n = Threads.atomic_add!(ndone, 1) + 1
                (n % 32 == 0 || n == P) && _log("    $out: $n/$P shards")
            end
            _mark(out)
            _log("  $out: $(Base.format_bytes(_dir_bytes(od))) on disk, $(Base.format_bytes(_free_bytes(work))) free on out_dir")
            for i in ins
                rm(joinpath(work, i); recursive = true, force = true)
            end
            if !tripped[] && startswith(out, "super_L1_")
                tripped[] = true
                pb = _partial_pool_bytes(od, P)                  # unpacked (sample,count) bytes here
                raw = pb * (final_total / covered[out]) / 1e9
                # V4.0 packs the pool block-FOR (~5-8x) and never expands it globally.
                _log("pool projection from $out: ~$(round(raw; digits = 1)) GB unpacked / " *
                     "~$(round(raw / 6; digits = 1)) GB packed (est.) at $final_total samples")
            end
        end

        # Fire every merge node whose inputs are all complete, repeatedly, so a finished
        # L1 immediately unblocks its L2 parent when the siblings were already done.
        function _drain_ready!()
            progressed = true
            while progressed
                progressed = false
                for (out, ins) in plan
                    _done(out) && continue
                    all(_done(i) for i in ins) || continue
                    _run_merge_node!(out, ins)
                    progressed = true
                end
            end
        end

        _drain_ready!()  # resume: collapse whatever the previous run left ready

        for w in 1:n_waves
            lo = (w - 1) * wave_size + 1
            hi = min(w * wave_size, length(sample_paths))
            first_global = n_prefix + lo
            name = "wave_$(w)"; d = joinpath(work, name)
            if !(resume && _done(name))
                rm(d; recursive = true, force = true); mkpath(d)
                _await_disk(work, margin, "wave $w")
                # tmp_dir holds the whole wave's records before folding (~620 GB for 50
                # GTEx samples at K=12 aa on oni's ~1.2 TB local SSD), so size this check
                # from the last wave's records, kept in a sidecar so a resume knows it too
                recf = joinpath(mark_dir, "wave_rec_bytes.n")
                last_rec = isfile(recf) ? parse(Int, read(recf, String)) : 0
                _await_disk(rec_root, round(Int, 1.1 * last_rec) + 20_000_000_000, "wave $w records (tmp_dir)")
                _log("wave $w/$n_waves: samples $(lo)..$(hi) (global $(first_global)..$(n_prefix + hi))")
                rec_bytes = _run_wave(sample_paths[lo:hi], first_global, d, joinpath(rec_root, name),
                          keep_shift, P, counter; count_attempts = count_attempts)
                write(joinpath(mark_dir, "wave_rec_bytes.n"), string(rec_bytes))
                rm(joinpath(rec_root, name); recursive = true, force = true)
                _mark(name)
            end
            _drain_ready!()
        end
        _drain_ready!()

        # reshard: split any shard whose combined size across every root exceeds
        # `reshard_leaf_bytes` further, by `reshard_extra_bits` additional k-mer prefix
        # bits, so assembly's concurrency window never has to hold several giant shards
        # at once (see ARCHITECTURE.md). `reshard_extra_bits` fixes the bit position every
        # split level consumes, so it's pinned the first time this stage runs against this
        # `work` dir and checked on every later (resumed) call, the same way a seed's
        # sample count is pinned via its own `.n` sidecar above.
        root_dirs = String[joinpath(work, r) for r in roots]
        ebfile = joinpath(mark_dir, "reshard_extra_bits.n")
        if isfile(ebfile)
            prior_eb = parse(Int, read(ebfile, String))
            prior_eb == reshard_extra_bits ||
                error("reshard_extra_bits=$reshard_extra_bits doesn't match this work dir's " *
                      "prior run ($prior_eb); use the same value to resume, or delete $work")
        else
            write(ebfile, string(reshard_extra_bits))
        end
        leaves = _reshard_leaves(root_dirs, P, keep_shift, work, reshard_leaf_bytes,
                                 reshard_extra_bits, _done, _mark)

        # assemble
        nk_sum = 0; np_sum = 0
        for lp in leaves, f in lp
            nk, np = _partial_sizes(f); nk_sum += nk; np_sum += np
        end
        _await_disk(dirname(out_path), round(Int, 2.5 * np_sum + 4.1 * nk_sum) + margin, "assembly")
        _log("assembling -> $out_path  ($final_total samples, $(length(leaves)) leaves)")
        _assemble_v4(leaves, kk, Ab, idx_prefix, final_total, out_path; nk_hint = nk_sum)
        rm(work; recursive = true, force = true)
        return out_path
    finally
        rm(rec_root; recursive = true, force = true)
    end
end
