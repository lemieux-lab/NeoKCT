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

const _STREAM_BITS_PER_AA = Int(bits_per_symbol(AAAlphabet()))  # 5

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

# Runs one wave: counts `paths` (global sample indices first_global .. first_global+len-1),
# streams records into per-shard files under `rec_dir`, then folds each shard into
# `out_dir/shard{p}.bin`. `counter(path)::Dict{UInt64,UInt32}`.
function _run_wave(paths::Vector{String}, first_global::Int, out_dir::String, rec_dir::String,
                   keep_shift::Int, P::Int, counter)
    mkpath(rec_dir)
    mask = UInt64(P - 1)
    flush_at = 1_000_000
    ios = IOStream[open(joinpath(rec_dir, "rec$(p).bin"), "w") for p in 0:(P - 1)]
    bufs = [_StreamRec[] for _ in 0:(P - 1)]
    try
        for (li, path) in enumerate(paths)
            slot = UInt32(first_global - 1 + li)
            ht = counter(path)::Dict{UInt64, UInt32}
            for (k, c) in ht
                p = _shard_of(k, keep_shift, mask) + 1
                push!(bufs[p], _StreamRec(k, slot, c))
                if length(bufs[p]) >= flush_at
                    write(ios[p], bufs[p]); empty!(bufs[p])
                end
            end
            ht = nothing
        end
        for p in 1:P
            !isempty(bufs[p]) && write(ios[p], bufs[p])
        end
    finally
        for io in ios; close(io); end
    end

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

function _assemble_v4(root_dirs::Vector{String}, P::Int, aak::Int, total_samples::Int, out_path::String)
    all_kmers = UInt64[]
    all_row_id = UInt32[]
    grp = _RowPool()

    for p in 0:(P - 1)
        shards = _PartialShard[_read_partial_shard(joinpath(d, "shard$(p).bin")) for d in root_dirs]
        m = _merge_partial_shards(shards)
        append!(all_kmers, m.kmers)
        @inbounds for c in eachindex(m.kmers)
            r = m.row_id[c]
            pairs = Tuple{UInt32, UInt32}[]
            for t in m.offsets[r]:(m.offsets[r + 1] - 1)
                push!(pairs, (m.pool_s[t], m.pool_c[t]))
            end
            push!(all_row_id, _intern_row!(grp, pairs))
        end
    end

    seqs = DeltaArray{UInt64, UInt32}(all_kmers, DEFAULT_CHECKPOINT_INTERVAL)  # already globally sorted
    kl = KmerLayer{aak, AAAlphabet, UInt64, UInt32}(seqs, _empty_kmer_idx(aak, AAAlphabet(), UInt64))
    scl = SparseCountsLayer(all_row_id, grp.offsets, grp.pool_s, grp.pool_c, Ref(Int64(total_samples)))
    kct = KCT(kl, scl)
    compute_index!(kct.kmer)
    write_kct(kct, out_path)
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

function _seed_row_pairs(scl::SparseCountsLayer, ::Nothing, i::Int, base::Int)
    r = scl.row_id[i]
    pairs = Tuple{UInt32, UInt32}[]
    @inbounds for t in scl.row_offsets[r]:(scl.row_offsets[r + 1] - 1)
        push!(pairs, (UInt32(base + scl.row_samples[t]), scl.row_counts[t]))
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
whole build.

# Keyword Arguments
- `K`: nucleotide k-mer length (amino-acid length is `K ÷ 3`). `chunks` is the read-chunk size
- `wave_size`: samples counted and folded per wave
- `shard_bits`: gives `2^shard_bits` prefix shards. `merge_fanin` partials are combined per
  hierarchical merge step
- `counter`: `counter(path)::Dict{UInt64,UInt32}` produces one sample's k-mer counts
  (default `jello_superthreaded_hash`)
- `seeds`: `Vector{Tuple{String,Int}}` of `(existing .kct, first_global_sample)` folded in as
  level-0 partials. Sample ranges must be disjoint and precede the waves
- `resume`: skips any wave or partial that already has a `DONE` marker
- `max_pool_gb`: aborts after the first super-partial if the projected row pool exceeds it
"""
function build_kct_streaming(sample_paths::Vector{String};
        K::Int = 30, chunks::Int = 500_000, wave_size::Int = 100, shard_bits::Int = 8,
        merge_fanin::Int = 5, max_pool_gb::Real = 400,
        tmp_dir::String = _STREAM_TMP_DEFAULT, out_dir::String,
        out_path::String = joinpath(out_dir, "neokct_v4.kct"),
        counter = (p -> jello_superthreaded_hash(p, K, chunks)),
        seeds::Vector{Tuple{String, Int}} = Tuple{String, Int}[],
        resume::Bool = true)

    aak = K ÷ 3
    P = 1 << shard_bits
    keep_shift = _STREAM_BITS_PER_AA * aak - shard_bits
    keep_shift >= 1 || error("shard_bits=$shard_bits too large for $(aak)-AA k-mers")
    merge_fanin >= 2 || error("merge_fanin must be >= 2")

    work = joinpath(out_dir, "streaming_work"); mkpath(work)
    rec_root = joinpath(tmp_dir, "neokct_stream_rec_$(getpid())"); mkpath(rec_root)
    _done(name) = isfile(joinpath(work, name, "DONE"))
    _mark(name) = touch(joinpath(work, name, "DONE"))
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
                _mark(name)
            end
            nS = open(io -> (read(io, Int64); read(io, Int64)), joinpath(d, "manifest"))
            covered[name] = nS; n_prefix += nS
            push!(seed_names, name)
            println("seed $si: $path -> $nS samples (global $(first_s)..$(n_prefix))")
        end

        # waves
        n_waves = cld(length(sample_paths), wave_size)
        wave_names = String[]
        total_samples = n_prefix
        for w in 1:n_waves
            lo = (w - 1) * wave_size + 1
            hi = min(w * wave_size, length(sample_paths))
            first_global = n_prefix + lo
            name = "wave_$(w)"; d = joinpath(work, name)
            if !(resume && _done(name))
                rm(d; recursive = true, force = true); mkpath(d)
                println("wave $w/$n_waves: samples $(lo)..$(hi) (global $(first_global)..$(n_prefix + hi))")
                _run_wave(sample_paths[lo:hi], first_global, d, joinpath(rec_root, name),
                          keep_shift, P, counter)
                rm(joinpath(rec_root, name); recursive = true, force = true)
                _mark(name)
            end
            covered[name] = hi - lo + 1
            total_samples = n_prefix + hi
            push!(wave_names, name)
        end

        # hierarchical merge
        plan, roots = _merge_plan(vcat(seed_names, wave_names), merge_fanin)
        checked = false
        for (out, ins) in plan
            covered[out] = sum(covered[i] for i in ins)
            od = joinpath(work, out)
            if !(resume && _done(out))
                rm(od; recursive = true, force = true); mkpath(od)
                println("merge $out <- [$(join(ins, ", "))]  ($(covered[out]) samples)")
                for p in 0:(P - 1)
                    shards = _PartialShard[_read_partial_shard(joinpath(work, i, "shard$(p).bin")) for i in ins]
                    _write_partial_shard(joinpath(od, "shard$(p).bin"), _merge_partial_shards(shards))
                end
                _mark(out)
                for i in ins; rm(joinpath(work, i); recursive = true, force = true); end
            end
            if !checked && startswith(out, "super_L1_")
                checked = true
                pb = _partial_pool_bytes(od, P)
                proj = pb * (total_samples / covered[out]) / 1e9
                println("row-pool tripwire: $(round(pb / 1e9; digits = 2)) GB over $(covered[out]) samples " *
                        "-> ~$(round(proj; digits = 1)) GB projected at $total_samples samples")
                proj > max_pool_gb && error("projected row pool ~$(round(proj; digits = 1)) GB exceeds " *
                    "max_pool_gb=$max_pool_gb; row diversity is higher than expected (see plan fallback)")
            end
        end

        # assemble
        println("assembling -> $out_path  ($total_samples samples)")
        _assemble_v4([joinpath(work, r) for r in roots], P, aak, total_samples, out_path)
        rm(work; recursive = true, force = true)
        return out_path
    finally
        rm(rec_root; recursive = true, force = true)
    end
end
