# NeoKCT architecture

Full design rationale for the codebase, kept separate from the source so the inline comments can stay short. If a comment in the code says "see ARCHITECTURE.md", the detail is here.

## Layered design

A `KCT{K, Ab, Counts, Biotype, C, D}` is a `KmerLayer` plus up to two optional layers, `Counts` and `Biotype`, each either a real layer type or `Nothing`. The compiler dispatches `getindex`, `add_biotypes`, etc. on the concrete type parameters, so there is no runtime layer check. This replaced three separate historical types (`NeoKCT`, `RichKCT`, `GenomicIndex`), which are now just particular `Counts`/`Biotype` combinations.

- `K`: k-mer length in the layer's own alphabet symbols (amino acids for `AAAlphabet`, nucleotides for `DNAAlphabet{2}`).
- `Ab`: the alphabet.
- `Counts`: `CountsLayer` (streaming-built, posting-list-backed) or `Nothing`.
- `Biotype`: `BiotypLayer` or `Nothing`.
- `C`, `D`: checkpoint and delta word types for the `KmerLayer`'s `DeltaArray`.

An earlier incremental design (`push!` / `collapse!` growing a CSR/`PackedArray`-backed `CountsLayer` one batch at a time) existed through V3.0 and was removed once the streaming builder made it obsolete for every real use case: rebuilding from scratch with `build_kct_streaming` is faster than even a single incremental `push!` in practice, and the old structure's per-k-mer index grew superlinearly with sample count (78% of a ~150 GB table at 140 GTEx samples). Its file formats (V1.2 through V3.0, and RichKCT V2.0) are no longer readable; check out the `V4.0-last-retrocompat` git tag to read one of those files.

## KmerLayer: the sorted k-mer store

`KmerLayer.seqs` is a `DeltaArray{C,D}` (from the `NArrays` package): the sorted distinct k-mers as periodic full checkpoints (type `C`, default `UInt64`) plus small deltas between consecutive entries (type `D`, default `UInt32`). Sequential iteration costs one addition per k-mer. Random access costs an addition per entry back to the nearest checkpoint, so it is O(`DEFAULT_CHECKPOINT_INTERVAL` = 256). A delta that would overflow `D` is transparently promoted to a checkpoint, so any two consecutive stored entries are never more than a `D` step apart no matter how the k-mers are distributed. This halves storage vs. a raw `Vector{UInt64}` in practice, since lexicographically adjacent 64-bit k-mer codes tend to be numerically close.

### Prefix search index

`KmerLayer.idx` partitions the sorted array by the leading `prefix_size` symbols (default `DEFAULT_IDX_PREFIX_SIZE = 5`) into `1 << ((K - prefix_size) * bits_per_symbol(Ab))` buckets. Each bucket stores `(decoded k-mer value at the bucket's start, index range)`, so `findfirst` computes the bucket directly from the query's leading bits, then calls `searchfirst` (a linear delta-walk in `NArrays`, not a binary search) seeded at the bucket's known start value and position, skipping the need to rewind to a `DeltaArray` checkpoint.

`prefix_size` trades index memory against per-lookup scan length: each unit shrinks the bucket table by `2^bits_per_symbol` (32x for 5-bit AA) but grows the average bucket (and so the scan). At K=10 AA over ~6e9 k-mers, prefix 4 gives 2^30 buckets (~26 GB, ~5 k-mers/bucket) vs. prefix 5's 2^25 (~0.8 GB, ~180 k-mers/bucket); 5 is the chosen default because the index stays sub-GB and the extra scan length is negligible next to the rest of a lookup. DNA tables use a different derivation (see `build_kct_streaming` below) since a 30-symbol DNA k-mer has a far larger symbol-count range than a 10-symbol AA one.

`compute_index!` walks the sorted array once, recording where each distinct prefix group starts. The very first group (position 1) is sealed specially rather than via the normal "previous key" bookkeeping: seeding the walk with a `key=0` placeholder would look like a real seal for whatever k-mer's decoded value happens to be exactly 0 (e.g. an all-first-symbol k-mer), so a `seeded` flag distinguishes "haven't started" from "last key really was 0".

## CountsLayer: block-compressed sparse-count posting lists

Built only by `build_kct_streaming` (see below); read-only afterward (`push!`, `collapse!`, `sort!`, `repack` all raise `_COUNTS_READONLY`).

### Why this shape, not a dedup pool

An earlier design (row-deduplicated: one pointer per k-mer into a pool of distinct sparse rows) was measured against real GTEx data before being replaced. On a 500-sample AA build, rows longer than 4 pairs held ~96% of all pairs but deduplicated at ~1.0x (essentially unique), while short rows (<=4 pairs) deduplicated well (~12x) but only ever held ~4% of the pairs. A pool plus a per-k-mer `row_id` pointer therefore cost about as many bytes as it saved, on top of needing a global dedup dictionary sized to the number of distinct rows, tens of billions at full cohort scale, which does not fit in memory during assembly and was the direct cause of an OOM. Dropping cross-k-mer dedup entirely removes that dictionary and turns assembly into a pure streaming transcode.

### On-disk / in-memory layout

- `block_size::Int32`: k-mers per block (`_V4_BLOCK_KMERS = 64`).
- `block_ptr::Vector{UInt64}`: 0-based byte offset of each block in `blob`.
- `blob::Vector{UInt8}`: the packed blocks, concatenated. Memory-mapped on load (`Mmap`), not read into RAM, since at full cohort it can be 1-2 TB and a query touches one block; falls back to a plain `read!` if the filesystem refuses `mmap`.
- `n_kmers`, `n_samples`: `RefValue{Int64}`, stored explicitly since row data is inline (there is no pool whose length would otherwise imply the k-mer count).

Each k-mer's row is its `(sample, count)` pairs for the samples it was observed in, ascending by sample, with absent samples implicit. A block packs `block_size` consecutive k-mers' rows with one shared bit width per field, computed as the minimum that fits every row in the block (frame-of-reference / PFOR style):

```
block layout in blob, starting at block_ptr[b]+1 (1-based):
  UInt8 nk | UInt8 lbits | UInt8 gbits | UInt8 cbits | varint n_exc
  MSB-first bitstream: nk row lengths (lbits each)
                     | T sample gaps  (gbits each), T = sum(row lengths)
                     | T counts       (cbits each)
  byte-aligned, then n_exc records: varint local_pair_index, UInt32LE value
```

`lbits == 0` means every row in the block has length 1 (an all-singleton block); `cbits == 0` means every count is 1. Sample indices are gap-coded (each stored as the delta to the previous sample in the same row, first entry raw), which packs tight when the sample list is grouped so a k-mer's presence clusters (e.g. GTEx ordered by tissue). Counts above 65535 would force every row in the block to a wide `cbits`; instead the packed value saturates to the sentinel `0xffff` and the true value goes in a short per-block exception list, keyed by position within the block's flat pair stream (`_V4_COUNT_INLINE_MAX`). Measured: ~0.001% of pairs are exceptions.

Reading one row (`_decode_row` / `_decode_block`) decodes the whole 64-row block it belongs to; there is no cheaper single-row path, since 64 rows is little enough to decode in full each time and this keeps the codec simple. A bulk pass (`_each_row`) decodes each block once in order rather than once per k-mer.

### Bit codec (`_BitWriter` / `_BitReader`)

MSB-first, byte-oriented: `_BitWriter` accumulates into a `UInt64` and emits whole bytes as they fill (never more than 7 live bits carried between calls, since every `_put!` call here is <=16 bits), `_BitReader` mirrors it. `_bitlen` (`Base.top_set_bit`) gives the minimum width for a block's fields. Varints are standard LEB128, used for `n_exc` and exception indices since both are usually 0 or small.

## BiotypLayer

Per-k-mer biotype (protein_coding, lncRNA, intergenic, ...) as an interned bitmask: `ids` is one `UInt16` pool index per k-mer, `pool` is the deduplicated set of bitmasks actually seen (bit `b-1` of a mask corresponds to `biotype_names[b]`), `pool[1]` is always `INTERGENIC_MASK` and is the default for any k-mer with no annotation. `build_genomic_index` builds a `KCT{K÷3, AAAlphabet, Nothing, BiotypLayer}` from an Ensembl transcript FASTA plus a GTF/GFF3 (biotype from the FASTA header's `transcript_biotype`, falling back to the annotation file; GFF3 transcript IDs are Ensembl-prefixed, `"transcript:ENST..."`, and are stripped before matching). `add_biotypes(kct, gidx)` left-joins a table's k-mers against a genomic index by an O(n+m) sorted-merge walk over both `DeltaArray`s (both are already sorted, so this is a single linear pass, no per-k-mer binary search), producing a new KCT that carries the matched biotype alongside `kct`'s counts layer unchanged (the join only touches `kmer`/`biotype`).

## `.kct` file format

Every file starts with a leading `Float64` version tag, read by `load_kct` to pick the right reader (`_load_kct(io, Val(version))`). The only version `write_kct` produces now is 4.0 (`KCT_VERSION`); earlier versions are retained only as retrocompat history in the note above, not as writable or readable formats in this codebase.

```
Header:      [Float64 version][Int64 K][Int64 Ab_name_len][UInt8... Ab_name][UInt8 layers_mask][Int64 n_kmers]
KmerLayer:   [Int64 sizeof(C)][Int64 sizeof(D)][Int64 cp_interval][Int64 n_cp][Int64 n_rci][C... cps][D... deltas][Int64... rci]
CountsLayer  (layers_mask bit 0): [Int64 n_samples][Int64 n_kmers][Int32 block_size][Int64 n_blocks][Int64 blob_len][UInt64... block_ptr][UInt8... blob]
BiotypLayer  (layers_mask bit 1): [Int64 n_names]([Int64 len][UInt8... name]...)[Int64 pool_len][UInt64... pool][UInt16... ids]
```

`Ab_name` is `string(Ab)` (e.g. `"DNAAlphabet{2}"`, keeping type parameters, unlike `.name.singletonname`), parsed back with `eval(Meta.parse(...))`. `_write_counts_header` is factored out from `_write_counts` so `build_kct_streaming`'s assembly can write the header, then stream a multi-hundred-GB blob in from a temp file, without ever holding it in one `Vector` (see below). A retrocompat reader still exists for the standalone GenomicIndex v1.0 `.gidx` format (`load_gidx`), which never touches a counts layer and so was unaffected by the CountsLayer rename/removal; it loads into the current `KCT` type.

## `build_kct_streaming` (`StreamBuild.jl`): disk-backed external merge sort

The only builder in this codebase; alphabet-generic: `translate=false` (default) builds `KCT{K, DNAAlphabet{2}}`, `translate=true` builds `KCT{K÷3, AAAlphabet}`; everything below the counter works on raw `UInt64` k-mer codes regardless, and only the prefix-shard bit shift and the final `KmerLayer` construction need to know which. `idx_prefix` (the `KmerLayer` prefix index size) defaults differently per alphabet when unset: AA keeps the tuned `DEFAULT_IDX_PREFIX_SIZE = 5`, DNA uses `kk - 12` so a 30-symbol k-mer's bucket table still stays sub-GB (prefix 18 -> 2^24 buckets, ~400 MB).

This is an external k-way merge sort (Knuth TAOCP Vol 3 Section 5.4; same family as KMC and DSK): samples are counted in waves, each wave's observations are bucketed into 256 (by default, `shard_bits`) k-mer-prefix shards on local scratch disk, folded into sorted per-shard partials, and partials are combined by a hierarchical k-way merge before a final assembly pass. Nothing dense is ever resident; peak memory is bounded by one wave's working set plus whatever is in flight in the current merge step, flat in cohort size, until the final assembly.

### Pipeline stages

1. **Count, per wave** (`_run_wave`): `wave_size` (default 100) samples counted one at a time (`counter(path)::Dict{UInt64,UInt32}`, default `jello_superthreaded_hash`), each sample's records (`_StreamRec{kmer,slot,count}`, `slot` a 1-based *global* sample index) scattered by k-mer prefix (`_shard_of`) into `P = 2^shard_bits` record files on `tmp_dir` (local disk, not NFS: this is the highest-churn I/O of the whole build). `_count_with_retry` retries a sample up to `count_attempts` times with linear backoff before giving up and naming it; this exists because a single transient failure (busy-NFS short read, samtools SIGPIPE at EOF, a FIFO reader deadlocked on a samtools that died before opening it) used to kill a multi-day build outright with no exception for the caller to catch. A per-sample count that returns zero k-mers should be treated as suspect (a stub/truncated input), not a legitimate empty sample; callers building their own `counter` should raise in that case so retry has something to catch.

2. **Fold** (`_build_wave_shard`): each shard's records are sorted by k-mer (`psort!`), equal-k-mer runs collapsed into `(sample, count)` rows sorted by sample, and identical rows *within this shard* deduplicated via `_RowPool` (a `Dict` from the row's flattened `Vector{UInt32}` key to a pool index). This dedup is cheap here (rows are narrow, at most `wave_size` samples wide) and keeps the on-disk partial smaller; it is deliberately not carried through to the final assembly (see CountsLayer above). The result is a `_PartialShard` (sorted k-mers, `row_id` per k-mer, CSR pool).

3. **Hierarchical merge, interleaved with the wave loop**: partials are combined `merge_fanin` (default 5) at a time via `_merge_plan` (a bottom-up plan of `super_L{level}_{group}` nodes with deterministic names, independent of *when* each runs) and `_merge_partial_shards` (a k-way merge across the group, same per-shard row dedup as the fold step). Critically, `_drain_ready!` fires every plan node whose inputs are already marked done, called after every wave (and once before the wave loop, for resume), not only after all waves finish. This is what bounds live partial storage on disk to O(`merge_fanin`) waves instead of O(all waves): the naive "merge everything only at the end" design was tried first and needs O(n_waves) partials resident simultaneously, which exhausted a shared NFS export mid-build in production before this fix.

   Markers live in a flat `work/_markers/` directory, not inside each partial's own folder: once merges are interleaved, a wave's folder is deleted (by the merge that consumes it) long before the whole build finishes, and a marker kept inside it would vanish with the folder, making `resume` recount an already-merged wave. A completed leaf's or node's `covered` (sample count) is likewise computed once, up front, as a pure function of the wave/seed list (`n_waves`, `wave_size`, and each seed's own recorded count survive in `_markers/*.n` files), so it never needs to re-derive itself from a folder that might already be gone.

4. **Assembly** (`_assemble_v4`): the terminal step, the only one where a real `KCT` / `CountsLayer` is constructed. For each of the `P` shards (processed in a bounded parallel window, `W = min(nthreads, 8)`, since NFS read latency benefits from concurrency but each in-flight merge holds several partials' worth of memory), `_merge_shard_rows` k-way-merges that shard across the final `<= merge_fanin` root partials **without any deduplication**, in flat CSR form (`kmers`, `off`, `flat_s`, `flat_c`, no `Vector` built per k-mer). The roots cover disjoint, ascending global-sample ranges and every stored row is already sample-sorted, so a k-mer's full row across roots is just its per-root pieces concatenated in root order, no re-sort needed. (An earlier version of this step built one small `Vector{UInt32}` per k-mer per root; at ~300M k-mers in a single dense shard that is several hundred million tiny allocations, which was the dominant cost once the old cross-root dedup dictionary had already been removed. Flat CSR replaces that with a handful of `append!`-grown arrays per shard.)

   The main thread consumes each shard's rows in order, buffering up to `_V4_BLOCK_KMERS` (64) at a time before calling `_pack_block!` and appending the packed bytes to a growing in-memory buffer, which is itself flushed to a `blobtmp` file (next to `out_path`) every 64 MiB. The blob is **never held whole in RAM**: at full cohort scale it can be 1-2 TB, and RAM instead holds only `all_kmers` (freed immediately after the `DeltaArray` is built), the block-offset array, and the current in-flight merge window. The final `.kct` is written by streaming that temp file's bytes straight into place after the header and k-mer block, then deleting it.

### Seeding and resume

`stream_kct_as_shards` turns an existing `.kct` into a level-0 partial with its samples renumbered from a given global offset, so a prior table can be folded in as a `seeds` entry without recounting those samples; seed and wave global sample ranges must be disjoint and seeds must precede the waves. `resume=true` (default) makes every stage check its flat marker before redoing work, so a killed build (or a deliberately pruned sample list, e.g. to drop a corrupt input file after the fact) picks back up rather than restarting.

## JelloFish: k-mer counting

Two counting paths coexist:

- `translate` + `k_merize` (older): builds a `Vector{Kmer}` per read and constructs a fresh `Kmer` object per position. Correct and general (handles arbitrary alphabets via the `Kmers`/`BioSequences` machinery) but allocates on every k-mer, which dominates runtime at the ~1e10 read-k-mers/sample scale a production build reaches. Still used by `GenomicIndexBuilder.jl` (a one-time reference build, not per-sample-at-scale) and `JellyfishDump.jl`'s translation step.
- The rolling counter (`count_kmers`, current default for all per-sample counting): keeps a 2-bit-per-base code in a `UInt64` accumulator, shifting in one base at a time (`_NT2BIT`, a 256-entry ASCII-to-2-bit lookup, `0xff` for anything not A/C/G/T/a/c/g/t, checked once per read via `_has_nonacgt` so a read with any ambiguous base is skipped whole rather than crashing the chunk). With `translate=false` the rolling code itself is the key; with `translate=true`, `_fold_codons` fields off 6-bit codons from the high end of the window and looks each up in `_CODON2AA5` (a 64-entry codon-to-5-bit-AA table built directly from `BioSequences.ncbi_trans_table[1]`, so it is bit-identical to `translate(k_merize(...))` on the same input) to build a `5*(K/3)`-bit AA code; hitting a stop codon anywhere in the window drops that k-mer (`_STOP_KEY` sentinel). Nothing is allocated in the per-base inner loop.

`chunk_stream` reads a bio file (`BioParser.jl`) in fixed-size chunks and spawns one counting task per chunk, gated by a `Base.Semaphore(max_inflight)` (default `nthreads()`): the reader blocks once that many tasks are outstanding, rather than spawning unboundedly ahead of the counters, which used to grow an unbounded task backlog under sustained high read-in rate. `jello_superthreaded_hash` folds the per-chunk `Dict`s together in a binary tree of `merge` tasks running concurrently with the counting, not as a separate serial pass at the end.

## JellyfishDump: external Jellyfish binary format

Parses the *external* `jellyfish count`/`jellyfish dump` binary format (unrelated to this project's own `JelloFish.jl`), for cases where only a pre-computed `.jf` dump is available (e.g. archived GTEx counts) and raw reads are not. Layout: a decimal byte-offset, a `{`-delimited JSON header (cmdline args recovered from it: `-m`/`--mer-len` for K, `--out-counter-len` for the per-record counter width), then fixed-size `(k-mer, count)` records. When a k-mer's `K*2` bits and the counter width together still fit one 64-bit word, both are packed into a single record instead of two fields, mirroring how `Kmers.jl` itself right-aligns a k-mer in a `UInt64` (data in the low bits, zero padding above); `_clean_kmer_bits` masks off that padding so the raw word is a valid Dict key and sorts/compares correctly as a plain `UInt64` downstream. Canonical dumps (`-C`) are rejected by default (`allow_canonical=true` to override), since a canonical k-mer's stored bit pattern could be either strand and translating it in one fixed frame would silently produce the wrong peptide for about half of them. `verify_jellyfish_dump` cross-checks a prefix of the parsed records against the `jellyfish` CLI's own `dump -c` text output, for validating this parser against an unfamiliar Jellyfish build before trusting a full production parse (the binary record layout was reverse-engineered, not read from a format spec).

## AAAlphabet

A 5-bit `BioSequences`/`Kmers` alphabet: an `AminoAcid`'s own byte value (0..21, up to `AA_Gap`) is used directly as the encoding, so a K-mer of amino acids packs into `K*5` bits. Beyond the standard `encode`/`decode` hooks, `Kmers` needs bit-addressing helpers for pulling one symbol out of a multi-word k-mer whose symbols do not align to 64-bit boundaries: `offset` (unused high-bit padding in the first word), `true_index` (logical symbol index to absolute bit position), `chunk` (which `.data` word a symbol lands in), and `extract_encoded_element` (handles a symbol that straddles two words).

## Benchmarking (`KCTBenchmarker.jl`)

`benchmark_kct` measures a table's component sizes (k-mer sequence store, block index, packed blob), an extrapolated pair count (sampling up to 2000 blocks rather than decoding all of them), bytes/k-mer, and query speed, appending an entry to a JSON history (`benchmark_data_v4.json`). It times `findfirst` over a sample of the table's own k-mers (drawn from `DeltaArray` checkpoints, so all guaranteed hits), which measures k-mer *location*, not count-vector decode. There are no plots; an earlier version of this file drew CairoMakie growth-over-samples SVGs for the incremental `CountsLayer`'s component vocabulary (`flat_cids`, `n_cids`, packed words, bitmap), which no longer applies now that the incremental layer is gone.
