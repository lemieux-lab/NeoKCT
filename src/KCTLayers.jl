using Kmers
using BioSequences
using BioSymbols
using ProgressMeter
using Base.Threads
using NArrays
import NArrays: repack
export repack
using Dates
using BitIntegers

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

# Number of trailing k-mer symbols the prefix index does NOT partition on. The lookup
# table has `1 << ((K - DEFAULT_IDX_PREFIX_SIZE) * bits_per_symbol)` buckets, so every +1
# here shrinks it by `2^bits_per_symbol` (32x for a 5-bit AA alphabet). findfirst then
# linearly scans the bucket it lands in (searchfirst is a linear delta-walk, not a binary
# search), so this trades index RAM for a longer per-lookup scan: at K=10/AAAlphabet,
# 4 -> 1<<30 buckets (~26 GB, ~5 k-mers/bucket at 6e9 k-mers) vs 5 -> 1<<25 (~0.8 GB,
# ~180 k-mers/bucket). 5 keeps the table sub-GB while the scan stays short enough not to
# matter next to the rest of a push!.
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

    # Walk the sorted k-mer array and record where each distinct prefix group starts,
    # along with the decoded value at that position (the same O(1)-amortised pass regular
    # iteration already does), so findfirst can seed searchfirst directly instead of
    # rewinding to a DeltaArray checkpoint. On a prefix transition at position i, seal the
    # range for the previous key and open a new one from i, findfirst uses these ranges to
    # skip most of the array.
    #
    # The very first group (i=1) is handled separately from later transitions: sealing it
    # with the initial start_val=zero(C) placeholder (as if key=0 had really been seen)
    # would be wrong whenever a real k-mer's own decoded value is exactly 0 (e.g. an
    # all-A/all-first-symbol k-mer, a real possibility with e.g. poly-A runs). That
    # zero would then look like a spurious match for that placeholder-seeded bucket.
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

CSR-style store of per-k-mer count vectors across all loaded samples.

Each k-mer maps to one or more count words held in a `PackedArray`.
The mapping is encoded as a flat list of count-word IDs (`flat_cids`) and a
parallel `n_cids` vector recording how many IDs each k-mer owns, so offsets
can be reconstructed on the fly without storing them explicitly.

# Fields
- `n_cids`: number of count-word IDs per k-mer, cumulative sum reconstructs offsets
- `flat_cids`: flat array of count-word IDs, one contiguous block per k-mer
- `counts`: packed count words, deduplicated across k-mers with identical count vectors
- `samples`: number of samples currently stored in this layer
- `shared_words`: word IDs referenced by more than one k-mer, cached, see `find_shared_words`
"""
# Mutable so collapse! / _merge_and_sort! / sort! can hand over a freshly built `flat_cids`
# (or deduped `counts`) by rebinding the field. That is an O(1) swap that lets the old buffer
# be collected, instead of `resize!` + `copyto!` into the existing one, which needs both the
# old and the new array (~125 GB each at production scale) resident at the same time.
mutable struct CountsLayer <: AbstractLayer
    n_cids::Vector{UInt16}  # number of cids per k-mer (cumsum to reconstruct offset vector)
    flat_cids::Vector{UInt32}
    counts::PackedArray{UInt32}
    samples::Base.RefValue{Int64}
    shared_words::Set{UInt32}
end

Base.length(cl::CountsLayer) = length(cl.n_cids)

CountsLayer(; word_size::Type{<:Unsigned}=UInt128) =
    CountsLayer(UInt16[], UInt32[], PackedArray{UInt32, word_size}(), Ref(0), Set{UInt32}())

# Convenience constructor for callers that don't already know shared_words (e.g. loading
# a KCT from disk): (re)derives it with one find_shared_words scan. Fine here since this
# only runs once at load time, not once per sample, push! must use the cached field
# instead (see the note on find_shared_words).
CountsLayer(n_cids::Vector{UInt16}, flat_cids::Vector{UInt32}, counts::PackedArray{UInt32},
            samples::Base.RefValue{Int64}) =
    CountsLayer(n_cids, flat_cids, counts, samples,
                find_shared_words(CountsLayer(n_cids, flat_cids, counts, samples, Set{UInt32}())))

function assemble_count_vector(cl::CountsLayer, i::Integer)
    # lo re-sums every preceding n_cids on each call, O(i) per call, so a loop reading
    # every k-mer (e.g. kct[1:end]) is O(n^2). Fine for spot lookups. Would need cached
    # offsets (see _kmer_offsets) if bulk reads become a bottleneck.
    lo = 1 + Int64(sum(@view cl.n_cids[1:i-1]))
    cids = @view cl.flat_cids[lo : lo + cl.n_cids[i] - 1]
    counts = UInt32[]
    for c in cids
        append!(counts, cl.counts[Int(c)])
    end

    # Trailing zeros (samples where this k-mer was absent) are stored implicitly to save space.
    # Materialise them only when reading, by padding to the full sample count.
    n_missing = cl.samples.x - length(counts)
    return n_missing > 0 ? vcat(counts, zeros(eltype(cl.counts.words), n_missing)) : counts
end

Base.getindex(cl::CountsLayer, i::Integer) = assemble_count_vector(cl, i)

function find_shared_words(cl::CountsLayer)

    # After collapse!, k-mers with identical count vectors point to the same PackedArray
    # word IDs. push! must not append to a shared word or it would corrupt every other
    # k-mer that references it. Identify all word IDs referenced more than once.
    #
    # This is an O(flat_cids) full-table scan: only used to (re)establish shared_words
    # from scratch (loading from disk, after a repack). CountsLayer caches the result on
    # its `shared_words` field otherwise. push! must never call this itself, since
    # collapse! is the only operation that creates sharing (see _shared_from_global_perm)
    # and push! runs once per sample, which would make this an O(samples * table size)
    # cost across a build.
    #
    # Use a count array indexed by word ID: O(words) space, much smaller than a Set over
    # all flat_cids when words >> unique k-mers (e.g. before collapse, 4 bytes per word
    # vs ~50 bytes per k-mer for the old Set{UInt32} approach).
    n_words = length(cl.counts.words)
    ref_counts = zeros(UInt8, n_words)  # UInt8: saturate at 2, we only need ≤ 2
    for cid in cl.flat_cids
        @inbounds ref_counts[cid] < 2 && (ref_counts[cid] += UInt8(1))
    end
    return Set{UInt32}(i for i in eachindex(ref_counts) if ref_counts[i] > 1)
end

# Derives the shared-word set as a side effect of collapse!'s own dedup pass, instead of
# a separate O(flat_cids) find_shared_words scan. global_perm[i] = the deduped word that
# original word i maps to. A deduped word is shared if EITHER:
#   - more than one original word maps to it this round (including two words from the
#     same k-mer's own chain, where appending to either would still corrupt the other, so
#     that must count as shared too, same as find_shared_words), or
#   - exactly one original word X maps to it, but X was already shared before this
#     collapse! (old_shared, from the previous round). Sharing is really a property of
#     how many flat_cids entries point at a word, and collapse! preserves that reference
#     count under the new id even when it doesn't dedupe that word against anything new
#     this round. Missing this case is what made the cached set silently shrink on every
#     later collapse!, since a word shared two rounds ago would just look "unmerged" (count
#     1) in this round's tally on its own.
function _shared_from_global_perm(global_perm::Vector{UInt32}, n_deduped::Int, old_shared::Set{UInt32})
    ref_counts = zeros(UInt8, n_deduped)
    for j in global_perm
        @inbounds ref_counts[j] < 2 && (ref_counts[j] += UInt8(1))
    end
    shared = Set{UInt32}(i for i in eachindex(ref_counts) if ref_counts[i] > 1)
    for x in old_shared
        push!(shared, global_perm[x])
    end
    return shared
end

function _kmer_offsets(n_cids::Vector{UInt16})::Vector{UInt64}
    offsets = Vector{UInt64}(undef, length(n_cids) + 1)
    offsets[1] = 1
    for i in eachindex(n_cids)
        offsets[i+1] = offsets[i] + n_cids[i]
    end
    return offsets
end

function collapse!(cl::CountsLayer)
    deduped, _, global_perms = permdedup(cl.counts)
    # Remap every cid to its deduped word id in place: each entry is read then written at
    # the same index and global_perms is read-only, so the threads touch disjoint slots.
    # Avoids the full `similar(cl.flat_cids)` copy the old version allocated (~125 GB at
    # production scale, held alongside the original).
    @threads for i in eachindex(cl.flat_cids)
        @inbounds cl.flat_cids[i] = UInt32(global_perms[cl.flat_cids[i]])
    end
    cl.shared_words = _shared_from_global_perm(global_perms, length(deduped.words), cl.shared_words)
    cl.counts = deduped
    printstyled("Collapse done\n", color=:green)
    return cl
end

# Rebuilds cl.counts with word type W2, preserving old word boundaries (a word is only
# split if its own values no longer fit under W2). Cheap, but a larger W2 can't shrink
# flat_cids since old words are never merged together.
function _repack_split(cl::CountsLayer, ::Type{W2}) where {W2<:Unsigned}
    new_counts, id_map = repack(cl.counts, W2)

    offsets = _kmer_offsets(cl.n_cids)
    new_flat = UInt32[]
    sizehint!(new_flat, length(cl.flat_cids))
    new_n_cids = similar(cl.n_cids)

    @showprogress "Rewriting Count IDs..." for i in eachindex(cl.n_cids)
        lo, hi = offsets[i], offsets[i+1] - 1
        start = length(new_flat)
        for cid in @view cl.flat_cids[lo:hi]
            append!(new_flat, id_map[cid])
        end
        new_n_cids[i] = UInt16(length(new_flat) - start)
    end

    # Splitting an overflowing word doesn't break sharing (every k-mer that referenced it
    # gets the same id_map[cid] sequence), so a previously-shared word is still shared
    # here, just under (possibly several) new IDs. Re-derive rather than remap id-by-id
    # (the 4-arg CountsLayer constructor below does this via find_shared_words). This
    # only runs once per repack call, not once per sample, so the O(table size) cost is fine.
    return CountsLayer(new_n_cids, new_flat, new_counts, Ref(cl.samples.x))
end

# Rebuilds cl.counts with word type W2 by repacking each k-mer's own values from
# scratch into as few W2 words as possible, ignoring old word boundaries. Unlike
# _repack_split, this lets a larger W2 shrink flat_cids. It also breaks any
# word-sharing from a prior collapse! (identical count vectors now live in separate
# fresh words), which is why repack() re-collapses by default when merge=true.
function _repack_merge(cl::CountsLayer, ::Type{W2}) where {W2<:Unsigned}
    new_counts = PackedArray{UInt32, W2}()
    new_flat = UInt32[]
    new_n_cids = similar(cl.n_cids)
    offsets = _kmer_offsets(cl.n_cids)

    @showprogress "Repacking Counts Words (merging)..." for i in eachindex(cl.n_cids)
        lo, hi = offsets[i], offsets[i+1] - 1
        # append! is amortized O(1) per element. reduce(vcat, ...) here would recopy
        # everything gathered so far at each step, O(m^2) in the k-mer's own word count m.
        vals = UInt32[]
        for c in @view cl.flat_cids[lo:hi]
            append!(vals, cl.counts[Int(c)])
        end
        start = length(new_flat)
        if !isempty(vals)
            wid = new_word!(new_counts)
            push!(new_flat, UInt32(wid))
            last_set = 0
            for v in vals
                new_wid, last_set = push_at!(new_counts, v, wid, last_set)
                if new_wid != wid
                    wid = new_wid
                    push!(new_flat, UInt32(wid))
                end
            end
        end
        new_n_cids[i] = UInt16(length(new_flat) - start)
    end

    # Every k-mer gets a brand-new, isolated word chain here, so nothing is shared yet
    # (see repack()'s docstring: this is exactly the sharing this mode breaks).
    return CountsLayer(new_n_cids, new_flat, new_counts, Ref(cl.samples.x), Set{UInt32}())
end

# Lets word size be swept without re-parsing samples. `merge=true` consolidates each
# k-mer's values into as few new words as possible (shows the real benefit of a larger
# W2), `merge=false` (default) only splits words that overflow under W2, which is
# cheaper but can't shrink flat_cids. `collapse` re-dedupes afterwards and defaults to
# matching `merge`, since merging discards prior word-sharing from collapse!.
function repack(cl::CountsLayer, ::Type{W2}; merge::Bool=false, collapse::Bool=merge) where {W2<:Unsigned}
    new_cl = merge ? _repack_merge(cl, W2) : _repack_split(cl, W2)
    return collapse ? collapse!(new_cl) : new_cl
end

# V4.0 counts layer (row-deduplicated + sparse). Defined here so it can appear in the KCT
# `Counts` type parameter below. The streaming builder that produces it is in StreamBuild.jl.
include("SparseCountsLayer.jl")

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
specialised `getindex` and `collapse!` methods with no runtime overhead.

# Type Parameters
- `K`: k-mer length in amino-acid symbols
- `Ab`: alphabet type (e.g. `AAAlphabet`)
- `Counts`: `CountsLayer` (V3.0) or `SparseCountsLayer` (V4.0) when count data is present,
  `Nothing` otherwise
- `Biotype`: `BiotypLayer` when biotype data is present, `Nothing` otherwise
- `C`: checkpoint word type forwarded to the underlying `KmerLayer`
- `D`: delta word type forwarded to the underlying `KmerLayer`

# Fields
- `kmer`: sorted k-mer store
- `counts`: per-sample count layer, or `nothing`
- `biotype`: per-k-mer biotype annotation layer, or `nothing`
"""
struct KCT{K, Ab <: Alphabet, Counts <: Union{CountsLayer, SparseCountsLayer, Nothing}, Biotype <: Union{BiotypLayer, Nothing}, C <: Unsigned, D <: Unsigned}
    kmer::KmerLayer{K, Ab, C, D}
    counts::Counts
    biotype::Biotype
end

KCT(kl::KmerLayer{K, Ab, C, D}) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, Nothing, Nothing, C, D}(kl, nothing, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, Nothing, C, D}(kl, cl, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, BiotypLayer, C, D}(kl, cl, bl)
KCT(kl::KmerLayer{K, Ab, C, D}, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, Nothing, BiotypLayer, C, D}(kl, nothing, bl)
KCT(kl::KmerLayer{K, Ab, C, D}, scl::SparseCountsLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, SparseCountsLayer, Nothing, C, D}(kl, scl, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, scl::SparseCountsLayer, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, SparseCountsLayer, BiotypLayer, C, D}(kl, scl, bl)

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

function Base.getindex(kct::KCT{K, Ab, CountsLayer, Nothing}, i::Integer) where {K, Ab}
    return kct.kmer[i] => kct.counts[i]
end

function Base.getindex(kct::KCT{K, Ab, CountsLayer, BiotypLayer}, i::Integer) where {K, Ab}
    return kct.kmer[i] => (; counts=kct.counts[i], biotype=biotype_mask(kct.biotype, i))
end

function Base.getindex(kct::KCT{K, Ab, Nothing, BiotypLayer}, i::Integer) where {K, Ab}
    return kct.kmer[i] => biotype_mask(kct.biotype, i)
end

# V4.0 (SparseCountsLayer) reads: same shape as the CountsLayer methods above.
function Base.getindex(kct::KCT{K, Ab, SparseCountsLayer, Nothing}, i::Integer) where {K, Ab}
    return kct.kmer[i] => kct.counts[i]
end

function Base.getindex(kct::KCT{K, Ab, SparseCountsLayer, BiotypLayer}, i::Integer) where {K, Ab}
    return kct.kmer[i] => (; counts=kct.counts[i], biotype=biotype_mask(kct.biotype, i))
end

Base.getindex(kct::KCT, i::UnitRange) = Tuple(kct[j] for j in i)
Base.getindex(kct::KCT, i::AbstractVector) = Tuple(kct[j] for j in i)

"""
    collapse!(kct) -> KCT

Deduplicate the count layer of a V3.0 (`CountsLayer`) KCT: k-mers whose count
vectors are byte-identical are pointed at a single shared `PackedArray` word
chain, and `flat_cids` is remapped in place. Returns a `KCT` sharing the same
`kmer` (and `biotype`) layer with the collapsed counts. The build drivers call
this once per batch, right after `push!`. Errors on a V4.0 table.
"""
function collapse!(kct::KCT{K, Ab, CountsLayer, Nothing}) where {K, Ab}
    return KCT(kct.kmer, collapse!(kct.counts))
end

function collapse!(kct::KCT{K, Ab, CountsLayer, BiotypLayer}) where {K, Ab}
    return KCT(kct.kmer, collapse!(kct.counts), kct.biotype)
end

# V4.0 tables are assembled whole by build_kct_streaming and never mutated afterwards.
const _V4_READONLY = "V4.0 (SparseCountsLayer) KCTs are read-only; build them with build_kct_streaming"
collapse!(::KCT{K, Ab, SparseCountsLayer}) where {K, Ab} = error(_V4_READONLY)
Base.sort!(::KCT{K, Ab, SparseCountsLayer}) where {K, Ab} = error(_V4_READONLY)
Base.push!(::KCT{K, Ab, SparseCountsLayer}, ::Any) where {K, Ab} = error(_V4_READONLY)
repack(::KCT{K, Ab, SparseCountsLayer}, ::Type{<:Unsigned}; kwargs...) where {K, Ab} = error(_V4_READONLY)

"""
    repack(kct, W2; merge=false, collapse=merge) -> KCT

Rebuild the whole table's counts under packed-word type `W2`, for example to
sweep word sizes on an already-built KCT without re-parsing the original samples.
`merge=false` only splits words that overflow under `W2` (cheap, cannot shrink
`flat_cids`). `merge=true` repacks each k-mer's values from scratch into as few
`W2` words as possible, which can shrink `flat_cids` but drops the word-sharing
from a prior `collapse!`, so `collapse` re-runs afterwards and defaults to
matching `merge`. Errors on a V4.0 table.
"""
repack(kct::KCT{K, Ab, CountsLayer, Nothing}, ::Type{W2}; merge::Bool=false, collapse::Bool=merge) where {K, Ab, W2<:Unsigned} =
    KCT(kct.kmer, repack(kct.counts, W2; merge, collapse))

repack(kct::KCT{K, Ab, CountsLayer, BiotypLayer}, ::Type{W2}; merge::Bool=false, collapse::Bool=merge) where {K, Ab, W2<:Unsigned} =
    KCT(kct.kmer, repack(kct.counts, W2; merge, collapse), kct.biotype)

## KCT Sample Building ##

function KCT{K, Ab}(sample_hashtable::Dict{UInt64, UInt32};
                    checkpoint_size::Type{C}=UInt64,
                    delta_size::Type{D}=UInt32,
                    word_size::Type{<:Unsigned}=UInt128) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned}
    kl = KmerLayer{K, Ab}(checkpoint_size=checkpoint_size, delta_size=delta_size)
    cl = CountsLayer(word_size=word_size)
    tmp_seqs = Vector{UInt64}(undef, length(sample_hashtable))
    sizehint!(cl.flat_cids, length(sample_hashtable))
    sizehint!(cl.n_cids, length(sample_hashtable))

    # Pack each count into its own fresh word (one k-mer per word at this stage).
    # k-mers arrive in hash-table order, so store sequences in a temporary array for sorting.
    @showprogress desc="Parsing Hash-Table into KCT..." for (j, (k_bits, count)) in enumerate(sample_hashtable)
        push!(cl.counts, UInt64(count))
        tmp_seqs[j] = k_bits
        push!(cl.flat_cids, UInt32(lastindex(cl.counts)))
        push!(cl.n_cids, UInt16(1))
    end

    # Sort k-mers and permute the CSR arrays to restore sorted order.
    perm = psortperm(tmp_seqs)
    new_flat = similar(cl.flat_cids); new_n_cids = similar(cl.n_cids)
    for (j, i) in enumerate(perm)
        new_flat[j] = cl.flat_cids[i]
        new_n_cids[j] = cl.n_cids[i]
    end
    copyto!(cl.flat_cids, new_flat); copyto!(cl.n_cids, new_n_cids)
    encode!(kl.seqs, tmp_seqs[perm])
    cl.samples.x = 1
    kct = KCT(kl, cl)
    compute_index!(kct.kmer)
    return kct
end

## KCT push! helpers ##

# One (k-mer, batch-slot, count) observation. Batched push! flattens the whole batch of
# sample hash-tables into a Vector of these and sorts by k-mer, so each distinct k-mer is
# located in the table once per batch (not once per sample) and the O(table-size) merge +
# index rebuild + collapse run once per batch instead of once per sample.
struct _KmerCountRec
    kmer::UInt64
    slot::UInt16  # 1-based index of the sample within the batch
    count::UInt32
end
Base.isless(a::_KmerCountRec, b::_KmerCountRec) = isless(a.kmer, b.kmer)

# Append a value stream (`missing_counts` implicit-zero backfills, then this batch's counts
# up to `last_nz`) onto word `wid`, recording every additional word spilled into via `rec`.
# `rec` is a callback so the two call sites can differ only in how the extension list is
# looked up (eagerly for a shared last word, lazily otherwise). See the callers.
@inline function _append_count_stream!(cl::CountsLayer, wid::Int, W::Type,
                                       missing_counts::Int, cvec::AbstractVector{UInt32},
                                       last_nz::Int, rec::F) where {F}
    @inbounds for s in 1:(missing_counts + last_nz)
        v = s <= missing_counts ? W(0) : W(cvec[s - missing_counts] & typemax(W))
        nw = push_to!(cl.counts, v, wid)
        if nw != wid
            rec(nw)
            wid = nw
        end
    end
    return wid
end

# Brand-new k-mer (absent from every prior sample): fresh isolated word chain, one implicit
# zero per prior sample backfilled, then this batch's counts up to the last non-zero slot
# (`last_nz`). Slots after that stay implicit trailing zeros.
function _push_new_kmer_counts_batch!(cl::CountsLayer, cvec::AbstractVector{UInt32},
                                      last_nz::Int, prior_samples::Int)
    W = eltype(cl.counts.words)
    wid = new_word!(cl.counts)
    chunk_ids = UInt32[wid]
    wid = _append_count_stream!(cl, wid, W, prior_samples, cvec, last_nz,
                                nw -> push!(chunk_ids, UInt32(nw)))
    return chunk_ids
end

# Existing k-mer: `missing_counts` implicit-zero backfills (samples between this k-mer's
# last stored count and now) followed by this batch's counts up to `last_nz`. Mirrors the
# single-count path, appending in place when the k-mer's last word is unshared, otherwise
# starting a fresh chain so a word shared via collapse! is never mutated.
function _push_existing_kmer_counts_batch!(cl::CountsLayer, ext_buf::Dict{Int, Vector{UInt32}},
                                           k_pos::Int, shared_words::Set{UInt32},
                                           cvec::AbstractVector{UInt32}, last_nz::Int,
                                           offsets::Vector{UInt64})
    W = eltype(cl.counts.words)
    cur_cids = @view cl.flat_cids[offsets[k_pos] : offsets[k_pos+1]-1]
    last_wid = Int(cur_cids[end])
    total_stored = sum(length(cl.counts[Int(c)]) for c in cur_cids)
    missing_counts = cl.samples.x - total_stored

    if last_wid in shared_words
        # Shared last word: appending in place would corrupt the other k-mers that
        # reference it. Start a fresh, unshared chain, its first word going in ext eagerly.
        ext = get!(ext_buf, k_pos, UInt32[])
        wid = new_word!(cl.counts)
        push!(ext, UInt32(wid))
        _append_count_stream!(cl, wid, W, missing_counts, cvec, last_nz, nw -> push!(ext, UInt32(nw)))
    else
        # Unshared last word: append directly, touching ext_buf only on an actual spill
        # (the common no-backfill, count-fits case allocates nothing).
        _append_count_stream!(cl, last_wid, W, missing_counts, cvec, last_nz,
                              nw -> push!(get!(ext_buf, k_pos, UInt32[]), UInt32(nw)))
    end
    return
end

function _merge_and_sort!(kl::KmerLayer, cl::CountsLayer, bl::Union{BiotypLayer, Nothing},
                           ext_buf::Dict{Int, Vector{UInt32}}, new_seqs::Vector{UInt64},
                           new_cids::Vector{Vector{UInt32}}, offsets::Vector{UInt64})
    n_existing = length(kl.seqs)
    n_new = length(new_seqs)
    n_total = n_existing + n_new

    # kl.seqs is already sorted and no new k-mer equals an existing one (findfirst would
    # have matched it in push!), so sort only the small new batch (O(m log m)) and merge
    # the two disjoint sorted runs in one O(n_total) pass. kl.seqs is streamed through its
    # DeltaArray iterator rather than `collect`ed into a flat Vector first. That copy was
    # ~8 bytes * n_total (~45 GB at production scale) of pure transient garbage.
    new_perm = psortperm(new_seqs)

    merged_seqs = Vector{UInt64}(undef, n_total)
    new_flat = UInt32[]
    sizehint!(new_flat, length(cl.flat_cids) + sum(length, values(ext_buf); init=0) +
                        sum(length, new_cids; init=0))
    new_n_cids = Vector{UInt16}(undef, n_total)
    new_ids = isnothing(bl) ? nothing : Vector{UInt16}(undef, n_total)

    # Walk both sorted runs in lockstep, rebuilding flat_cids and n_cids in one pass.
    old_state = iterate(kl.seqs)  # nothing iff kl.seqs is empty
    oi = 1  # 1-based position of the old k-mer old_state currently holds
    ni = 1  # pointer into new_perm (i.e. into new_seqs, in sorted order)
    for j in 1:n_total
        cids_start = length(new_flat)
        old_val = isnothing(old_state) ? typemax(UInt64) : old_state[1]
        if ni > n_new || old_val < new_seqs[new_perm[ni]]

            # Existing k-mer: copy its original cid block, then any extensions from ext_buf.
            merged_seqs[j] = old_val
            append!(new_flat, @view cl.flat_cids[offsets[oi] : offsets[oi+1]-1])
            haskey(ext_buf, oi) && append!(new_flat, ext_buf[oi])
            isnothing(new_ids) || (new_ids[j] = bl.ids[oi])
            oi += 1
            old_state = isnothing(old_state) ? nothing : iterate(kl.seqs, old_state[2])
        else

            # Brand-new k-mer: use the word-ID list built by _push_new_kmer_counts_batch!.
            i = new_perm[ni]
            merged_seqs[j] = new_seqs[i]
            append!(new_flat, new_cids[i])
            isnothing(new_ids) || (new_ids[j] = UInt16(1))  # intergenic for new k-mers
            ni += 1
        end
        new_n_cids[j] = UInt16(length(new_flat) - cids_start)
    end

    if !isnothing(bl)
        resize!(bl.ids, n_total)
        copyto!(bl.ids, new_ids)
    end
    encode!(kl.seqs, merged_seqs)
    merged_seqs = UInt64[]  # free ~8 bytes * n_total before handing over the big cid array

    # O(1) handoff (CountsLayer is mutable): the old flat_cids / n_cids become garbage
    # instead of being kept alive for a resize! + copyto! into them.
    cl.flat_cids = new_flat
    cl.n_cids = new_n_cids
    return
end

## KCT push! and sort! ##

# Add one sample. Thin wrapper over the batched push! (a batch of one).
Base.push!(kct::KCT{K, Ab, CountsLayer, Biotype},
           sample_hashtable::AbstractDict{UInt64, UInt32}) where {K, Ab<:Alphabet, Biotype} =
    push!(kct, [sample_hashtable])

# Add a batch of samples in one shot. Locates each distinct k-mer once, appends all of the
# batch's counts for it, then runs the O(table-size) sorted merge + index rebuild a single
# time for the whole batch. Callers still collapse! afterwards (also once per batch).
# Returns the new total sample count. Mutates kct in place and does not return the KCT.
function Base.push!(kct::KCT{K, Ab, CountsLayer, Biotype},
                    batch::AbstractVector{<:AbstractDict{UInt64, UInt32}}) where {K, Ab<:Alphabet, Biotype}
    cl = kct.counts
    B = length(batch)
    B == 0 && return cl.samples.x

    # shared_words is cached on cl (collapse! is the only thing that creates sharing and it
    # isn't called here). Snapshot the CSR offsets before the loop. Those do change as we
    # append, but only into ext_buf / new_cids, never into cl.flat_cids itself, so the
    # snapshot stays valid for the whole loop and for _merge_and_sort!.
    shared_words = cl.shared_words
    offsets = _kmer_offsets(cl.n_cids)
    ext_buf = Dict{Int, Vector{UInt32}}()  # new word-ID extensions for existing k-mers
    new_seqs = UInt64[]  # bit-encodings of brand-new k-mers
    new_cids = Vector{UInt32}[]  # their word-ID lists

    # Flatten the batch into (k-mer, slot, count) records, sorted by k-mer so every
    # observation of a given k-mer across the batch's samples is contiguous.
    total = sum(length, batch; init=0)
    recs = Vector{_KmerCountRec}(undef, total)
    p = 0
    for slot in 1:B, (k, c) in batch[slot]
        p += 1
        @inbounds recs[p] = _KmerCountRec(k, UInt16(slot), c)
    end
    psort!(recs)

    prior_samples = cl.samples.x
    cbuf = Vector{UInt32}(undef, B)  # reused scratch: this k-mer's counts across the batch
    @showprogress desc="Adding Samples $(prior_samples+1)-$(prior_samples+B) to Table..." for r in 1:total
        (r > 1 && recs[r].kmer == recs[r-1].kmer) && continue  # only act on a run's first record
        k = recs[r].kmer
        fill!(cbuf, zero(UInt32))
        e = r
        @inbounds while e <= total && recs[e].kmer == k
            cbuf[recs[e].slot] = recs[e].count
            e += 1
        end
        last_nz = findlast(!iszero, cbuf)
        last_nz === nothing && continue  # k-mer with no non-zero count (shouldn't happen)

        k_pos = findfirst(kct.kmer, Kmer{Ab, K, 1}(Kmers.unsafe, (k,)))
        if k_pos == 0
            push!(new_seqs, k)
            push!(new_cids, _push_new_kmer_counts_batch!(cl, cbuf, last_nz, prior_samples))
        else
            _push_existing_kmer_counts_batch!(cl, ext_buf, k_pos, shared_words, cbuf, last_nz, offsets)
        end
    end

    recs = _KmerCountRec[]  # release (~16 bytes * total) before the merge allocates

    _merge_and_sort!(kct.kmer, cl, kct.biotype, ext_buf, new_seqs, new_cids, offsets)
    # The k-mer set (and therefore every prefix-index range) only changes when new k-mers
    # were added, so skip the full-table index scan otherwise.
    isempty(new_seqs) || compute_index!(kct.kmer)
    cl.samples.x += B
    return cl.samples.x
end

function Base.sort!(kct::KCT{K, Ab, CountsLayer, Biotype}) where {K, Ab, Biotype}
    isempty(kct.kmer.seqs) && return kct
    seqs = collect(kct.kmer.seqs)
    perm = psortperm(seqs)
    encode!(kct.kmer.seqs, seqs[perm])

    # Permute flat_cids and n_cids to match the new k-mer order.
    cl = kct.counts
    offsets = _kmer_offsets(cl.n_cids)
    new_flat = UInt32[]; sizehint!(new_flat, length(cl.flat_cids))
    new_n_cids = UInt16[]; sizehint!(new_n_cids, length(cl.n_cids))
    for i in perm
        lo, hi = offsets[i], offsets[i+1] - 1
        append!(new_flat, @view cl.flat_cids[lo:hi])
        push!(new_n_cids, cl.n_cids[i])
    end
    cl.flat_cids = new_flat  # O(1) handoff (mutable struct), see _merge_and_sort!
    cl.n_cids = new_n_cids

    if !isnothing(kct.biotype)
        new_ids = kct.biotype.ids[perm]
        copyto!(kct.biotype.ids, new_ids)
    end
    return kct
end

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

Works with either counts layer. The body only walks `kct.kmer.seqs` and forwards
`kct.counts` untouched, so a `CountsLayer` (V3.0) or `SparseCountsLayer` (V4.0)
table both pass through.
"""
function add_biotypes(kct::KCT{K, Ab, Counts, Nothing, C, D},
                      gidx::KCT{K, Ab, Nothing, BiotypLayer}) where {K, Ab, Counts<:Union{CountsLayer, SparseCountsLayer}, C<:Unsigned, D<:Unsigned}
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

"""
    build_kct(samples, K=30, chunks=500_000; kwargs...) -> KCT

Build a V3.0 (`CountsLayer`) KCT incrementally from a list of sample FASTQ paths.
The first sample seeds the table. The rest are counted in groups of `batch`, and
each group is merged into the table and `collapse!`d in one shot, so the
O(table-size) merge and dedup run once per batch instead of once per sample.
`K` is the nucleotide k-mer length (amino-acid length is `K ÷ 3`), `chunks` is
the read-chunk size handed to the parallel counter.

Checkpoints and benchmarks fire on batch boundaries. `save_at_samples` writes a
`.kct` at the first boundary at or past each listed sample count, into
`save_path`. A benchmark runs whenever the running total has advanced by
`benchmark_every` samples since the last one, writing plots and a JSON entry
under `save_path * "benchmarks/"`. `word_size`, `checkpoint_size` and
`delta_size` set the packed-count and DeltaArray word widths.
`full_pointer_walkthrough` adds a full `summarysize` pass to each benchmark.
"""
function build_kct(samples::AbstractVector{String}, K::Int=30, chunks::Int=500_000;
                   word_size::Type{<:Unsigned}=UInt128, checkpoint_size::Type{<:Unsigned}=UInt64,
                   delta_size::Type{<:Unsigned}=UInt32, save_at_samples::AbstractVector{Int}=Int[],
                   save_path::String="", batch::Int=10, benchmark_every::Int=5,
                   full_pointer_walkthrough::Bool=false)
    kct = KCT{K÷3, AAAlphabet}(jello_superthreaded_hash(popfirst!(samples), K, chunks),
                                checkpoint_size=checkpoint_size,
                                delta_size=delta_size,
                                word_size=word_size)
    kct = collapse!(kct)
    isempty(samples) && return kct
    mkpath(save_path * "benchmarks/")

    _ckpt(n) = write_kct(kct, save_path *
        "$(day(now()))_$(month(now()))_$(year(now()))_[$(n)_samples]_neokct.kct")

    pending = Dict{UInt64, UInt32}[]
    n = 1; last_bench = 1
    for (idx, sample) in enumerate(samples)
        push!(pending, jello_superthreaded_hash(sample, K, chunks))
        (length(pending) < batch && idx != length(samples)) && continue

        prev_n = n
        push!(kct, pending)
        kct = collapse!(kct)
        n += length(pending)
        empty!(pending)
        GC.gc()

        any(s -> prev_n < s <= n, save_at_samples) && _ckpt(n)
        if n - last_bench >= benchmark_every || idx == length(samples)
            benchmark_kct(kct, save_path * "benchmarks/",
                          full_pointer_walkthrough=full_pointer_walkthrough)
            last_bench = n
        end
    end
    return kct
end

include("KCTLoader.jl")

# From-scratch streaming builder: k-way merge of per-sample sorted k-mer streams into a
# V4.0 (SparseCountsLayer) table. Needs load_kct / write_kct from KCTLoader.jl above.
include("StreamBuild.jl")
