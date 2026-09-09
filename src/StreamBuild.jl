## Streaming k-way merge builder -> .kct V4.0 (SparseCountsLayer) ##
#
# From-scratch build for large sample sets that the incremental push!/collapse! path cannot
# reach. Samples are counted in waves. Each wave's k-mers are bucketed by prefix into shard
# record files, sorted, and folded into a prefix-sharded "partial" table (sorted k-mers plus
# a per-shard deduplicated pool of sparse (sample, count) rows). Partials are combined by a
# hierarchical k-way merge and the survivors are assembled into one V4.0 .kct.
#
# Nothing dense is ever resident: every stage streams one k-mer (one short row) at a time.
# Peak memory is a wave's k-mer hash-table + one shard, flat across the whole build.

# The streaming builder is alphabet-generic: everything downstream of the counter works on
# raw UInt64 k-mer codes. Only the prefix shard shift and the final KmerLayer construction
# need the alphabet, and both are derived from `translate` (see build_kct_streaming): DNA
# gives KCT{K, DNAAlphabet{2}}, translate gives KCT{K÷3, AAAlphabet}.
const _STREAM_BITS_PER_AA = 5   # bits_per_symbol(AAAlphabet())
const _STREAM_BITS_PER_NT = 2   # bits_per_symbol(DNAAlphabet{2}())

# println + flush: stdout is block-buffered when redirected to a file, so without the flush
# a multi-hour wave shows nothing in the log until the buffer fills or the process exits.
# Locked so the shard-parallel merge/assembly loops don't interleave half-lines.
const _LOG_LOCK = ReentrantLock()
_log(args...) = lock(_LOG_LOCK) do
    println("[", Dates.format(now(), "HH:MM:SS"), "] ", args...); flush(stdout)
end

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

function _write_partial_shard(path::String, ps::_PartialShard)
    open(path, "w") do io
        write(io, Int64(length(ps.kmers)))
        write(io, ps.kmers)
        write(io, ps.row_id)
        write(io, Int64(length(ps.offsets) - 1))  # n_rows
        write(io, ps.offsets)
        write(io, Int64(length(ps.pool_s)))  # n_pairs
        write(io, ps.pool_s)
        write(io, ps.pool_c)
    end
end

function _read_partial_shard(path::String)::_PartialShard
    open(path, "r") do io
        nk = read(io, Int64)
        kmers = Vector{UInt64}(undef, nk); read!(io, kmers)
        row_id = Vector{UInt32}(undef, nk); read!(io, row_id)
        nr = read(io, Int64)
        offsets = Vector{UInt64}(undef, nr + 1); read!(io, offsets)
        np = read(io, Int64)
        pool_s = Vector{UInt32}(undef, np); read!(io, pool_s)
        pool_c = Vector{UInt32}(undef, np); read!(io, pool_c)
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

# `counter(path)` reconstructs a CRAM through samtools + a FIFO, and a single transient
# hiccup there (short read off a busy NFS export, samtools SIGPIPE at EOF, FASTQ desync)
# used to kill a multi-day build. Retry the sample a few times before giving up, and name it
# when we do.
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

    _log("  folding $P shard record files into partials")
    for p in 0:(P - 1)
        rp = joinpath(rec_dir, "rec$(p).bin")
        _build_wave_shard(rp, joinpath(out_dir, "shard$(p).bin"))
        rm(rp; force = true)
    end
    open(joinpath(out_dir, "manifest"), "w") do io
        write(io, Int64(first_global)); write(io, Int64(length(paths)))
    end
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

## Assembly of the final .kct ##

const _V4_BLOCK_KMERS = 64  # k-mers per packed block in the SparseCountsLayer

# Merge shard `q` of every root partial into per-k-mer rows, WITHOUT deduplication, in flat
# CSR form: row i is `flat_s[off[i]:off[i+1]-1]` / `flat_c[...]`. The roots cover disjoint,
# ascending global-sample ranges and every stored row is already sample-sorted, so a k-mer's
# full row is its pieces concatenated in root order: no _RowPool, no Dict, no re-sort, and
# no Vector-per-k-mer (a shard here can hold ~300M k-mers -- that allocation storm was the
# bottleneck once the dedup dict was gone). `_merge_partial_shards` still does per-k-mer
# dedup at the wave/super level, where rows are narrow and it keeps the partials small.
function _merge_shard_rows(root_dirs::Vector{String}, q::Int)
    shards = _PartialShard[_read_partial_shard(joinpath(d, "shard$(q).bin")) for d in root_dirs]
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

# Transcode the merged root partials into one V4.0 table: for every k-mer in global sorted
# order, take its (sample, count) row and pack it inline into a block. No cross-k-mer dedup.
# The packed blob is streamed to a temp file next to `out_path` and copied into the .kct at
# the end, so RAM holds only the k-mer array, the block index and one merge window -- at
# full cohort the blob itself is ~1-2 TB. Per-shard merges run in a bounded parallel window;
# the pack loop is serial because block order = k-mer order.
function _assemble_v4(root_dirs::Vector{String}, P::Int, kk::Int, Ab::Type, idx_prefix::Int,
                     total_samples::Int, out_path::String)
    B = _V4_BLOCK_KMERS
    blobtmp = string(out_path, ".blob.tmp")
    rm(blobtmp; force = true)

    all_kmers = UInt64[]
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

        p = 0
        while p < P
            chunk = p:min(p + W - 1, P - 1)
            tasks = [Threads.@spawn _merge_shard_rows(root_dirs, q) for q in chunk]
            for t in tasks
                ks, roff, fs, fc = fetch(t)
                append!(all_kmers, ks)
                @inbounds for i in eachindex(ks)
                    a = roff[i]; b = roff[i + 1] - 1
                    npairs += (b - a + 1)
                    append!(pend_s, @view fs[a:b])
                    append!(pend_c, @view fc[a:b])
                    push!(pend_off, length(pend_s) + 1)
                    length(pend_off) - 1 == B && flush_block!()
                end
            end
            p += W
            _log("  assemble: $(min(p, P))/$P shards packed, $(length(all_kmers)) k-mers, " *
                 "$(Base.format_bytes(blob_len)) blob")
        end
        flush_block!()
        isempty(blobbuf) || write(bio, blobbuf)
    end

    n_kmers = length(all_kmers)
    seqs = DeltaArray{UInt64, UInt32}(all_kmers, DEFAULT_CHECKPOINT_INTERVAL)  # already globally sorted
    all_kmers = UInt64[]  # ~8*n_kmers freed before the write
    kl = KmerLayer{kk, Ab, UInt64, UInt32}(seqs, _empty_kmer_idx(kk, Ab(), UInt64; prefix_size = idx_prefix))
    scl = SparseCountsLayer(B, block_ptr, UInt8[], n_kmers, total_samples)  # blob stays on disk
    kct = KCT(kl, scl)
    compute_index!(kct.kmer; prefix_size = idx_prefix)
    _log("  assembled: $n_kmers k-mers, $npairs pairs, blob $(Base.format_bytes(blob_len)), " *
         "block_ptr $(Base.format_bytes(8 * length(block_ptr)))")

    open(out_path, "w") do io
        write(io, 4.0)
        _write_header_and_kmers(io, kct)
        _write_sparse_counts_header(io, total_samples, n_kmers, B, length(block_ptr), blob_len, block_ptr)
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

## Seed adapter: an existing .kct -> a level-0 partial ##

# Sparse (global_slot, count) pairs for k-mer i of a V3.0 CountsLayer, using precomputed
# CSR offsets (assemble_count_vector is O(i) per call, unusable in a full-table loop).
function _seed_row_pairs(cl::CountsLayer, offsets::Vector{UInt64}, i::Int, base::Int)
    vals = UInt32[]
    for c in @view cl.flat_cids[offsets[i]:(offsets[i + 1] - 1)]
        append!(vals, cl.counts[Int(c)])
    end
    pairs = Tuple{UInt32, UInt32}[]
    @inbounds for (s, v) in enumerate(vals)
        v != 0 && push!(pairs, (UInt32(base + s), v))
    end
    return pairs
end

# V4.0 source: decode k-mer i's inline row. Decodes its whole block, so a full seed walk is
# O(block_size) per k-mer; acceptable for the occasional seed fold.
function _seed_row_pairs(scl::SparseCountsLayer, ::Nothing, i::Int, base::Int)
    s, c = _decode_row(scl, i)
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
    n_src = counts isa CountsLayer ? Int(counts.samples.x) : Int(counts.n_samples.x)
    v3_offsets = counts isa CountsLayer ? _kmer_offsets(counts.n_cids) : nothing
    base = first_global_sample - 1

    kmers = [UInt64[] for _ in 0:(P - 1)]
    rids = [UInt32[] for _ in 0:(P - 1)]
    rps = [_RowPool() for _ in 0:(P - 1)]
    mask = UInt64(P - 1)

    it = iterate(kct.kmer.seqs); i = 1
    while it !== nothing
        kbits, st = it
        p = _shard_of(kbits, keep_shift, mask) + 1
        pairs = _seed_row_pairs(counts, v3_offsets, i, base)
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

## Pool-size tripwire ##

function _partial_pool_bytes(dir::String, P::Int)
    total = 0
    for p in 0:(P - 1)
        open(joinpath(dir, "shard$(p).bin"), "r") do io
            nk = read(io, Int64)
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

Build a V4.0 (`SparseCountsLayer`) `.kct` from `sample_paths` by a k-way streaming merge.
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
- `wave_size`: samples counted and folded per wave
- `shard_bits`: gives `2^shard_bits` prefix shards. `merge_fanin` partials are combined per
  hierarchical merge step
- `counter`: `counter(path)::Dict{UInt64,UInt32}` produces one sample's k-mer counts
  (default `jello_superthreaded_hash`, honouring `translate`)
- `seeds`: `Vector{Tuple{String,Int}}` of `(existing .kct, first_global_sample)` folded in as
  level-0 partials. Sample ranges must be disjoint and precede the waves
- `resume`: skips any wave or partial that already has a `DONE` marker
- `count_attempts`: retries per sample when `counter` throws (transient samtools/FIFO/NFS
  failures), with a linear backoff, before aborting and naming the sample
- `max_pool_gb`: accepted for CLI compatibility, no longer used (the build logs a size
  projection after the first super-partial but never aborts)
"""
function build_kct_streaming(sample_paths::Vector{String};
        K::Int = 30, chunks::Int = 500_000, translate::Bool = false,
        wave_size::Int = 100, shard_bits::Int = 8,
        merge_fanin::Int = 5, max_pool_gb::Real = 400, idx_prefix::Int = -1,
        count_attempts::Int = 4,
        tmp_dir::String = _STREAM_TMP_DEFAULT, out_dir::String,
        out_path::String = joinpath(out_dir, "neokct_v4.kct"),
        counter = (p -> jello_superthreaded_hash(p, K, chunks; translate = translate)),
        seeds::Vector{Tuple{String, Int}} = Tuple{String, Int}[],
        resume::Bool = true)

    Ab = translate ? AAAlphabet : DNAAlphabet{2}
    bps = translate ? _STREAM_BITS_PER_AA : _STREAM_BITS_PER_NT
    kk = translate ? K ÷ 3 : K
    # `idx_prefix` < 0 means derive it: AA keeps the tuned prefix of 5 (1<<25 buckets), DNA
    # uses kk-12 so the bucket table stays sub-GB (k=30 -> prefix 18 -> 1<<24, ~400 MB).
    idx_prefix < 0 && (idx_prefix = translate ? DEFAULT_IDX_PREFIX_SIZE : max(0, kk - 12))
    P = 1 << shard_bits
    keep_shift = bps * kk - shard_bits
    keep_shift >= 1 || error("shard_bits=$shard_bits too large for a $(kk)-symbol k-mer")
    merge_fanin >= 2 || error("merge_fanin must be >= 2")

    work = joinpath(out_dir, "streaming_work"); mkpath(work)
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
            rm(od; recursive = true, force = true); mkpath(od)
            _log("  merge $out <- [$(join(ins, ", "))]  ($(covered[out]) samples), $(Threads.nthreads()) threads")
            in_dirs = String[joinpath(work, i) for i in ins]
            ndone = Threads.Atomic{Int}(0)
            Threads.@threads for p in 0:(P - 1)
                _write_partial_shard(joinpath(od, "shard$(p).bin"), _merge_one_shard(in_dirs, p))
                n = Threads.atomic_add!(ndone, 1) + 1
                (n % 32 == 0 || n == P) && _log("    $out: $n/$P shards")
            end
            _mark(out)
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
                _log("wave $w/$n_waves: samples $(lo)..$(hi) (global $(first_global)..$(n_prefix + hi))")
                _run_wave(sample_paths[lo:hi], first_global, d, joinpath(rec_root, name),
                          keep_shift, P, counter; count_attempts = count_attempts)
                rm(joinpath(rec_root, name); recursive = true, force = true)
                _mark(name)
            end
            _drain_ready!()
        end
        _drain_ready!()

        # assemble
        _log("assembling -> $out_path  ($final_total samples)")
        _assemble_v4([joinpath(work, r) for r in roots], P, kk, Ab, idx_prefix, final_total, out_path)
        rm(work; recursive = true, force = true)
        return out_path
    finally
        rm(rec_root; recursive = true, force = true)
    end
end
