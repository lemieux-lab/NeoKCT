using Kmers
using BioSequences
using BioSymbols
using ProgressMeter
using Base.Threads
using NArrays
using Dates

abstract type AbstractLayer end

## K-mer Layer ##

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
- `idx`: prefix-length reference paired with a per-prefix `UnitRange` lookup table
"""
struct KmerLayer{K, Ab <: Alphabet, C <: Unsigned, D <: Unsigned} <: AbstractLayer
    seqs::DeltaArray{C, D}
    idx::Pair{Base.RefValue{Int64}, Vector{UnitRange{Int64}}}
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

KmerLayer{K, Ab}(; checkpoint_size::Type{C}=UInt64,
                   delta_size::Type{D}=UInt32) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KmerLayer{K, Ab, C, D}(DeltaArray{checkpoint_size, delta_size}(DEFAULT_CHECKPOINT_INTERVAL),
                            Ref(20) => fill(0:-1, 4^15))

function Base.findfirst(kl::KmerLayer{K, Ab}, key::Kmer{Ab, K}) where {K, Ab<:Alphabet}
    key_bits = key.data[1]
    idx_key = (key_bits >> (idx_prefix_size(kl) * bits_per_symbol(Ab()))) + 1
    r = kl.idx[2][idx_key]
    isempty(r) && return 0
    return searchfirst(kl.seqs, key_bits, r.start, r.stop)
end

Base.findfirst(kl::KmerLayer{K, Ab}, key::UInt64) where {K, Ab<:Alphabet} =
    findfirst(kl, Kmer{Ab, K, 1}(Kmers.unsafe, (key,)))

function compute_index!(kl::KmerLayer{K, Ab}; prefix_size::Int64=4) where {K, Ab<:Alphabet}
    start = 1
    last_key = 0x0000000000000000
    n = length(kl.seqs)
    prefix_shift = prefix_size * bits_per_symbol(Ab())

    # Walk the sorted k-mer array and record where each distinct prefix group starts.
    # On a prefix transition at position i, seal the range for the previous key and
    # open a new one from i; findfirst uses these ranges to skip most of the array.
    @showprogress "Computing Binary Search Index..." for (i, val) in enumerate(kl.seqs)
        key = val >> prefix_shift
        if key > last_key
            kl.idx[2][last_key+1] = start:i
            start = i
            last_key = key
        end
    end
    kl.idx[2][last_key+1] = start:n
    kl.idx[1].x = prefix_size
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
- `n_cids`: number of count-word IDs per k-mer; cumulative sum reconstructs offsets
- `flat_cids`: flat array of count-word IDs, one contiguous block per k-mer
- `counts`: packed count words, deduplicated across k-mers with identical count vectors
- `samples`: number of samples currently stored in this layer
"""
struct CountsLayer <: AbstractLayer
    n_cids::Vector{UInt16}  # number of cids per k-mer (cumsum to reconstruct offset vector)
    flat_cids::Vector{UInt32}
    counts::PackedArray{UInt32}
    samples::Base.RefValue{Int64}
end

Base.length(cl::CountsLayer) = length(cl.n_cids)

CountsLayer(; word_size::Type{<:Unsigned}=UInt128) =
    CountsLayer(UInt16[], UInt32[], PackedArray{UInt32, word_size}(), Ref(0))

function assemble_count_vector(cl::CountsLayer, i::Integer)
    lo = 1 + Int64(sum(@view cl.n_cids[1:i-1]))
    cids = @view cl.flat_cids[lo : lo + cl.n_cids[i] - 1]
    counts = reduce(vcat, [cl.counts[Int(c)] for c in cids])

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
    new_flat = similar(cl.flat_cids)
    @threads for i in eachindex(cl.flat_cids)
        new_flat[i] = UInt32(global_perms[cl.flat_cids[i]])
    end
    printstyled("Collapse done\n", color=:green)
    return CountsLayer(cl.n_cids, new_flat, deduped, cl.samples)
end

## Biotype Layer ##

const INTERGENIC_MASK = UInt64(1)  # bit 0 reserved for intergenic

"""
    BiotypLayer <: AbstractLayer

Per-k-mer biotype annotation stored as interned bitmasks.

Each bit position corresponds to one biotype in `biotype_names`: bit `b-1` set
means the k-mer originates from `biotype_names[b]`. Bitmasks are deduplicated
into a `pool` so `ids` stores only a `UInt16` index per k-mer. `pool[1]` is
always `INTERGENIC_MASK` (bit 0 only); k-mers with no annotation default to it.

# Fields
- `ids`: one pool index per k-mer (1-based into `pool`)
- `pool`: deduplicated bitmasks; `pool[1] == INTERGENIC_MASK`
- `biotype_names`: human-readable name per bit position (`biotype_names[b]` ↔ bit `b-1`)
"""
struct BiotypLayer <: AbstractLayer
    ids::Vector{UInt16}    # one per k-mer, index into pool (1-based)
    pool::Vector{UInt64}   # deduplicated biotype bitmasks; pool[1] always == INTERGENIC_MASK
    biotype_names::Vector{String}  # biotype_names[b] = name for bit b-1
end

# All-intergenic layer; biotype_names[1] must be "intergenic".
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

# Insert mask into pool if absent; return its 1-based UInt16 index.
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
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, Nothing, C, D}(kl, cl, nothing)
KCT(kl::KmerLayer{K, Ab, C, D}, cl::CountsLayer, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, CountsLayer, BiotypLayer, C, D}(kl, cl, bl)
KCT(kl::KmerLayer{K, Ab, C, D}, bl::BiotypLayer) where {K, Ab<:Alphabet, C<:Unsigned, D<:Unsigned} =
    KCT{K, Ab, Nothing, BiotypLayer, C, D}(kl, nothing, bl)

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

Base.getindex(kct::KCT, i::UnitRange) = Tuple(kct[j] for j in i)
Base.getindex(kct::KCT, i::AbstractVector) = Tuple(kct[j] for j in i)

function collapse!(kct::KCT{K, Ab, CountsLayer, Nothing}) where {K, Ab}
    return KCT(kct.kmer, collapse!(kct.counts))
end

function collapse!(kct::KCT{K, Ab, CountsLayer, BiotypLayer}) where {K, Ab}
    return KCT(kct.kmer, collapse!(kct.counts), kct.biotype)
end

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
        new_flat[j]   = cl.flat_cids[i]
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

function _push_new_kmer_counts!(cl::CountsLayer, count::UInt32)
    W = eltype(cl.counts.words)

    # Allocate a fresh word chain isolated from all existing k-mers.
    wid = new_word!(cl.counts)
    chunk_ids = UInt32[wid]

    # Backfill one zero per prior sample: this k-mer was absent from all of them.
    # When a word fills up, push! spills to a new word; track each new word ID,
    # as every word referenced by a k-mer must appear in its flat_cids entry.
    @inbounds for _ in 1:cl.samples.x
        push!(cl.counts, W(0), wid)
        if lastindex(cl.counts) != wid
            wid = lastindex(cl.counts)
            push!(chunk_ids, UInt32(wid))
        end
    end

    # Append the actual count for the current sample.
    push!(cl.counts, W(count & typemax(W)), wid)
    lastindex(cl.counts) != wid && push!(chunk_ids, UInt32(lastindex(cl.counts)))
    return chunk_ids
end

function _push_existing_kmer_counts!(cl::CountsLayer, ext_buf::Dict{Int, Vector{UInt32}},
                                      k_pos::Int, shared_words::Set{UInt32}, count::UInt32,
                                      offsets::Vector{UInt64})
    W = eltype(cl.counts.words)
    cur_cids = @view cl.flat_cids[offsets[k_pos] : offsets[k_pos+1]-1]
    last_wid = Int(cur_cids[end])

    # Total count slots stored across all words for this k-mer so far.
    total_stored = sum(length(cl.counts[Int(c)]) for c in cur_cids)

    # Trailing zeros for samples where the k-mer was absent are stored implicitly.
    # Now that a new count is being added they are no longer trailing, so they must
    # be materialised before the actual count can be appended.
    missing_counts = cl.samples.x - total_stored

    if last_wid in shared_words

        # The last word is shared with at least one other k-mer: appending in place
        # would corrupt that k-mer's data. Start a fresh, unshared word chain instead.
        # We always produce at least one new word here, so eagerly create the ext entry.
        ext = get!(ext_buf, k_pos, UInt32[])
        new_wid = new_word!(cl.counts)
        push!(ext, UInt32(new_wid))
        for _ in 1:missing_counts
            push!(cl.counts, W(0), new_wid)
            next_wid = lastindex(cl.counts)
            if next_wid != new_wid
                push!(ext, UInt32(next_wid)); new_wid = next_wid
            end
        end
        push!(cl.counts, W(count & typemax(W)), new_wid)
        word_id = lastindex(cl.counts)
        word_id != new_wid && push!(ext, UInt32(word_id))
    else

        # Last word is unshared: safe to append directly to it.
        # Only touch ext_buf on actual overflow or backfill overflow — the common path
        # (missing_counts == 0 and count fits in the current word) allocates nothing.
        for _ in 1:missing_counts
            push!(cl.counts, W(0), last_wid)
            new_wid = lastindex(cl.counts)
            if last_wid != new_wid
                push!(get!(ext_buf, k_pos, UInt32[]), UInt32(new_wid))
                last_wid = new_wid
            end
        end
        push!(cl.counts, W(count & typemax(W)), last_wid)
        word_id = lastindex(cl.counts)
        last_wid != word_id && push!(get!(ext_buf, k_pos, UInt32[]), UInt32(word_id))
    end
end

function _merge_and_sort!(kl::KmerLayer, cl::CountsLayer, bl::Union{BiotypLayer, Nothing},
                           ext_buf::Dict{Int, Vector{UInt32}}, new_seqs::Vector{UInt64},
                           new_cids::Vector{Vector{UInt32}}, offsets::Vector{UInt64})
    old_seqs = collect(kl.seqs)
    n_existing = length(old_seqs)

    # Merge existing and brand-new k-mer sequences, then sort to restore the invariant.
    all_seqs = vcat(old_seqs, new_seqs)
    perm = psortperm(all_seqs)

    new_flat = UInt32[]
    sizehint!(new_flat, length(cl.flat_cids) + sum(length, values(ext_buf); init=0) +
                        sum(length, new_cids; init=0))
    new_n_cids = UInt16[]
    sizehint!(new_n_cids, length(all_seqs))
    new_ids = isnothing(bl) ? nothing : Vector{UInt16}(undef, length(all_seqs))

    # Walk k-mers in sorted order, rebuilding flat_cids and n_cids in one pass.
    for (j, i) in enumerate(perm)
        cids_start = length(new_flat)
        if i <= n_existing

            # Existing k-mer: copy its original cid block, then any extensions from ext_buf.
            append!(new_flat, @view cl.flat_cids[offsets[i] : offsets[i+1]-1])
            haskey(ext_buf, i) && append!(new_flat, ext_buf[i])
            !isnothing(new_ids) && (new_ids[j] = bl.ids[i])
        else

            # Brand-new k-mer: use the word-ID list built by _push_new_kmer_counts!.
            append!(new_flat, new_cids[i - n_existing])
            !isnothing(new_ids) && (new_ids[j] = UInt16(1))  # intergenic for new k-mers
        end
        push!(new_n_cids, UInt16(length(new_flat) - cids_start))
    end

    if !isnothing(bl)
        resize!(bl.ids, length(all_seqs))
        copyto!(bl.ids, new_ids)
    end
    encode!(kl.seqs, all_seqs[perm])
    resize!(cl.flat_cids, length(new_flat)); copyto!(cl.flat_cids, new_flat)
    resize!(cl.n_cids, length(new_n_cids)); copyto!(cl.n_cids, new_n_cids)
end

## KCT push! and sort! ##

function Base.push!(kct::KCT{K, Ab, CountsLayer, Biotype},
                    sample_hashtable::Dict{UInt64, UInt32}) where {K, Ab<:Alphabet, Biotype}
    cl = kct.counts

    # Snapshot shared words and CSR offsets before the loop: both change as we append.
    shared_words = find_shared_words(cl)
    offsets = _kmer_offsets(cl.n_cids)
    ext_buf = Dict{Int, Vector{UInt32}}()  # new word-ID extensions for existing k-mers
    new_seqs = UInt64[]                    # bit-encodings of brand-new k-mers
    new_cids = Vector{UInt32}[]            # their word-ID lists

    @showprogress desc="Adding Sample $(cl.samples.x + 1) to Table..." for (k_bits, count) in sample_hashtable
        k_pos = findfirst(kct.kmer, Kmer{Ab, K, 1}(Kmers.unsafe, (k_bits,)))
        if k_pos == 0
            push!(new_seqs, k_bits)
            push!(new_cids, _push_new_kmer_counts!(cl, count))
        else
            _push_existing_kmer_counts!(cl, ext_buf, k_pos, shared_words, count, offsets)
        end
    end

    # Merge new k-mers into the sorted sequence array and rebuild the CSR structure.
    _merge_and_sort!(kct.kmer, cl, kct.biotype, ext_buf, new_seqs, new_cids, offsets)
    compute_index!(kct.kmer)
    cl.samples.x += 1
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
    resize!(cl.flat_cids, length(new_flat)); copyto!(cl.flat_cids, new_flat)
    resize!(cl.n_cids, length(new_n_cids)); copyto!(cl.n_cids, new_n_cids)

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
                                  Ref(20) => fill(0:-1, 4^15))
    kct = KCT(kl, BiotypLayer(ids, pool, biotype_names))
    compute_index!(kct.kmer)
    return kct
end

# Left-join kct k-mers against gidx via an O(n+m) sorted merge walk.
# Matched k-mers get the genomic biotype; unmatched k-mers stay intergenic.
function add_biotypes(kct::KCT{K, Ab, CountsLayer, Nothing, C, D},
                      gidx::KCT{K, Ab, Nothing, BiotypLayer}) where {K, Ab, C<:Unsigned, D<:Unsigned}
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
include("GenomicIndexBuilder.jl")
include("KCTBenchmarker.jl")

function build_kct(samples::AbstractVector{String}, K::Int=30, chunks::Int=500_000;
                   word_size::Type{<:Unsigned}=UInt128, checkpoint_size::Type{<:Unsigned}=UInt64,
                   delta_size::Type{<:Unsigned}=UInt32, save_at_samples::AbstractVector{Int}=Int[],
                   save_path::String="", collapse_every::Int=1, benchmark_every::Int=5,
                   full_pointer_walkthrough::Bool=false)
    kct = KCT{K÷3, AAAlphabet}(jello_superthreaded_hash(popfirst!(samples), K, chunks),
                                checkpoint_size=checkpoint_size,
                                delta_size=delta_size,
                                word_size=word_size)
    length(samples) != 0 && mkpath(save_path * "benchmarks/")
    for (i, sample) in enumerate(samples)
        push!(kct, jello_superthreaded_hash(sample, K, chunks))
        kct = ((i + 1) % collapse_every == 0) ? collapse!(kct) : kct
        (i + 1) in save_at_samples && write_kct(kct, save_path *
            "$(day(now()))_$(month(now()))_$(year(now()))_[$(i+1)_samples]_neokct.kct")
        ((i + 1) % benchmark_every == 0) && benchmark_kct(kct, save_path * "benchmarks/",
            full_pointer_walkthrough=full_pointer_walkthrough)
    end
    return kct
end

include("KCTLoader.jl")

# temporary path to quickly test on all samples
samples = readlines(open("/u/jacquinn/phd_stuff/data/tcga_fastqs.paths", "r"));