## Peptide queries ##

# Peptide-level lookups against an AA CountsLayer table (manuscript Methods 2.7). A peptide
# shorter than K resolves to a contiguous run of the sorted k-mer array: the first residue
# sits in the k-mer code's high bits (same layout as JelloFish's _fold_codons), so integer
# and lexicographic order coincide and every K-mer starting with a given prefix is adjacent.
# A peptide longer than K has no single matching k-mer, so its abundance is bounded above by
# the minimum count over its overlapping K-mer tiles.

# 5-bit AAAlphabet code of a peptide, first residue in the high bits, matching _fold_codons'
# layout so a translated read window and a peptide string hash to the same key.
function _aa_code(pep::AbstractString)
    code = UInt64(0)
    for c in pep
        e = BioSequences.encode(AAAlphabet(), AminoAcid(c))
        isnothing(e) && throw(ArgumentError("residue $c in $pep has no AAAlphabet code"))
        code = (code << 5) | UInt64(e)
    end
    return code
end

# (decoded value at position i, iterate state) such that iterate(a, st) resumes right after
# i, same walk as getindex but also handing back the state a sequential scan needs. Used
# once per prefix-range query to seed the walk at the bucket's start instead of re-decoding
# every position from its nearest checkpoint.
function _state_at(a::DeltaArray{C}, i::Int) where {C}
    inter = a.checkpoint_interval
    k = (i - 1) ÷ inter
    cp = a.regular_cp_idx[k + 1]
    v = a.checkpoints[cp]
    @inbounds for j in (k * inter + 2):i
        a.deltas[j] == 0 ? (cp += 1; v = a.checkpoints[cp]) : (v += a.deltas[j])
    end
    return v, (i + 1, cp, v)
end

"""
    searchsorted(kl::KmerLayer, lo_key::UInt64, hi_key::UInt64) -> UnitRange{Int}

Positions of every k-mer with `lo_key <= code <= hi_key`. Both keys must fall in the same
prefix-index bucket: true whenever the varying low bits of the range fit inside the
`idx_prefix_size(kl)` trailing symbols the index does *not* partition on, i.e. for a peptide
prefix of at least `K - idx_prefix_size(kl)` residues (always true for the callers below).
Throws rather than silently reading past a bucket boundary if that doesn't hold.

O(bucket size) in the worst case: a linear walk seeded at the bucket's start, same cost
class as `searchfirst`. Buckets average ~12,000 k-mers on the production GTEx table, cheap
next to the block decode a hit then does. A future version could binary-search within the
bucket's own range instead, if this walk ever shows up in a profile.
"""
function Base.searchsorted(kl::KmerLayer{K, Ab}, lo_key::UInt64, hi_key::UInt64) where {K, Ab<:Alphabet}
    shift = idx_prefix_size(kl) * bits_per_symbol(Ab())
    lo_key >> shift == hi_key >> shift || throw(ArgumentError("range spans several prefix buckets"))
    _, r = kl.idx[2][(lo_key >> shift) + 1]
    isempty(r) && return 1:0
    v, st = _state_at(kl.seqs, r.start)
    i = r.start
    first_hit = 0
    while true
        v > hi_key && break
        v >= lo_key && first_hit == 0 && (first_hit = i)
        i == r.stop && (i += 1; break)
        v, st = iterate(kl.seqs, st)
        i += 1
    end
    return first_hit == 0 ? (1:0) : (first_hit:(i - 1))
end

# All k-mers whose first L residues equal pep (L <= K). L == K degenerates to an exact
# match: free == 0, lo == hi, a single-k-mer (or empty) range.
function Base.searchsorted(kl::KmerLayer{K, AAAlphabet}, pep::AbstractString) where {K}
    L = length(pep)
    L > K && throw(ArgumentError("peptide longer than K = $K, tile it instead"))
    free = 5 * (K - L)
    lo = _aa_code(pep) << free
    return searchsorted(kl, lo, lo | ((UInt64(1) << free) - 1))
end

# Per-sample sum of the count rows in r, decoding each touched block once rather than once
# per row (assemble_count_vector's O(1)-row cost would mean re-decoding the same shared
# block many times over a wide prefix range).
function _sum_rows(cl::CountsLayer, r::UnitRange{Int})
    v = zeros(UInt32, Int(cl.n_samples.x))
    isempty(r) && return v
    B = Int(cl.block_size)
    for b in ((first(r) - 1) ÷ B):((last(r) - 1) ÷ B)
        for (li, (s, c)) in enumerate(_decode_block(cl, b))
            (b * B + li) in r || continue
            @inbounds for t in eachindex(s)
                v[s[t]] += c[t]
            end
        end
    end
    return v
end

#= Suffix anchoring, the fix for the C-terminal blind spot (PAPER_TODO.md A5 side note).
The prefix query needs the read to continue 3(K - L) stop-free nucleotides past the peptide,
so a peptide whose next codon is a stop is invisible: RLATAPSEK has a stop after it at 240 of
its 246 genomic locations and the prefix query saw 6 reads where BamQuery saw 518 over the C5
small cohort. The K-mers that end with the peptide aren't a contiguous range (the free
residues are in the high bits), but there are only 20^(K-L) of them, each an exact lookup:
20 for a 9-mer, 400 for an 8-mer. No new index, no table rebuild. =#
const _AA20 = "ACDEFGHIKLMNPQRSTVWY"  # what translate() can emit, stops are dropped at build time
const _MAX_SUFFIX_PAD = 3  # 20^3 = 8,000 lookups, beyond that suffix anchoring is skipped

# Per-sample sum of the rows of every K-mer X * pep, X over all 20^(K - L) residue strings.
function _suffix_counts(kct::KCT{K, AAAlphabet, CountsLayer}, pep::AbstractString) where {K}
    L = length(pep)
    pad = K - L
    pep_code = _aa_code(pep)
    aa_codes = [UInt64(BioSequences.encode(AAAlphabet(), AminoAcid(c))) for c in _AA20]
    v = zeros(UInt32, Int(kct.counts.n_samples.x))
    for m in 0:(20^pad - 1)
        x = UInt64(0)
        r = m
        for _ in 1:pad
            x = (x << 5) | aa_codes[r % 20 + 1]
            r ÷= 20
        end
        i = findfirst(kct.kmer, (x << (5 * L)) | pep_code)
        i == 0 || (v .+= kct.counts[i])
    end
    return v
end

# Was this table built with cut-and-pad windows (pad=true)? Pad and stop codes never come out
# of plain translation, and in a padded table they sit in a large share of k-mers (every read
# end, every stop), so ~4k evenly spaced k-mers settle it without a format flag. Cached per
# table object, the answer can't change for a loaded table.
const _PADDED_CACHE = IdDict{Any, Bool}()
const _PADDED_LOCK = ReentrantLock()
function _is_padded(kct::KCT{K, AAAlphabet}) where {K}
    lock(_PADDED_LOCK) do
        get!(_PADDED_CACHE, kct) do
            n = length(kct.kmer.seqs)
            n == 0 && return false
            for i in unique(round.(Int, range(1, n; length=min(n, 4096))))
                v = kct.kmer.seqs[i]
                for j in 0:(K - 1)
                    c = (v >> (5j)) & 0x1f
                    (c == _PAD_AA5 || c == _STOP_AA5) && return true
                end
            end
            return false
        end
    end
end

"""
    peptide_counts(kct, pep; il_ambiguous=false, anchor=:auto) -> Vector{UInt32}

Per-sample abundance of `pep` against a `CountsLayer` table. `length(pep) <= K` sums the
counts over the peptide's whole prefix range: the number of read windows that start with a
coding sequence of `pep` and continue stop-free for `3(K - length(pep))` more nucleotides,
mirroring BamQuery's requirement that a read span the whole MCS. `length(pep) == K` is the
k-mer's own row. `length(pep) > K` takes the elementwise minimum over the
`length(pep) - K + 1` overlapping K-mer tiles, an upper bound on the number of read windows
carrying the whole peptide. `il_ambiguous=true` sums over every I/L variant of `pep`, for a
peptide whose I/L assignment came from a database lookup rather than the spectrum (leucine
and isoleucine are isobaric, indistinguishable by mass alone).

On a cut-and-pad table (built with `pad=true`) the prefix range alone is an exact count of
read occurrences for any `length(pep) <= K`: every occurrence contributes exactly one window
starting at it, cut and padded if it hit a stop, an N or the read's end. `anchor=:auto`
(default) uses `:prefix` there and `:both` on any other table.

`anchor` controls how a peptide shorter than K is anchored inside its K-mers. `:prefix`
is the original query above. `:both` also counts the read windows that *end*
with the peptide, every `X * pep` with `X` one of the `20^(K - length(pep))` residue strings
(one exact lookup each), and returns the per-sample maximum of the two. That recovers a
peptide sitting at the C-terminus of its ORF, whose prefix range is empty because the
codon right after it is a stop, and reads where the peptide sits at the read's own 3' edge.
The max is still a lower bound on the reads carrying the peptide and never counts one twice.
Suffix anchoring is skipped when `K - length(pep) > $(_MAX_SUFFIX_PAD)` (more than
$(20^_MAX_SUFFIX_PAD) lookups), which never happens for an 8-11-mer against K = 10.

Remaining blind spots: a peptide that is a whole ORF by itself (stop on both sides, no
room for either extension) still returns 0 for `length(pep) < K`; matching is exact, so a
read with a sequencing error inside the window is not counted, same as BamQuery.
"""
function peptide_counts(kct::KCT{K, AAAlphabet, CountsLayer}, pep::AbstractString;
                        il_ambiguous::Bool=false, anchor::Symbol=:auto) where {K}
    anchor in (:auto, :prefix, :both) || throw(ArgumentError("anchor must be :auto, :prefix or :both, got :$anchor"))
    anchor == :auto && (anchor = _is_padded(kct) ? :prefix : :both)
    il_ambiguous && return sum(peptide_counts(kct, v; anchor=anchor) for v in _il_variants(pep))
    L = length(pep)
    if L <= K
        v = _sum_rows(kct.counts, searchsorted(kct.kmer, pep))
        (anchor == :prefix || L == K || K - L > _MAX_SUFFIX_PAD) && return v
        return max.(v, _suffix_counts(kct, pep))
    end
    v = fill(typemax(UInt32), Int(kct.counts.n_samples.x))
    for s in 1:(L - K + 1)
        i = findfirst(kct.kmer, _aa_code(pep[s:s + K - 1]))
        i == 0 && return zeros(UInt32, length(v))
        v .= min.(v, kct.counts[i])
    end
    return v
end

# Every I/L variant of pep: one bit per I/L position, 2^n_IL strings (n_IL is small for a
# real MHC-I ligand, so this never gets close to enumerating anything large).
function _il_variants(pep::AbstractString)
    pos = findall(c -> c == 'I' || c == 'L', pep)
    out = String[]
    for m in 0:(2^length(pos) - 1)
        c = collect(pep)
        for (b, p) in enumerate(pos)
            c[p] = (m >> (b - 1)) & 1 == 1 ? 'I' : 'L'
        end
        push!(out, String(c))
    end
    return out
end

"""
    peptide_matrix(kct, peps; il_ambiguous=false, anchor=:auto) -> Matrix{UInt32}

`peptide_counts` for every peptide in `peps`, one row per peptide, one column per sample.
Peptides are independent, so the loop runs across threads. Logs progress through `_Prog`
(Progress.jl), since a real peptide set against the full GTEx table is exactly the kind of
run PAPER_TODO.md's conventions want instrumented.
"""
function peptide_matrix(kct::KCT{K, AAAlphabet, CountsLayer}, peps::Vector{String};
                        il_ambiguous::Bool=false, anchor::Symbol=:auto) where {K}
    M = zeros(UInt32, length(peps), Int(kct.counts.n_samples.x))
    p = _Prog(length(peps), "peptide_matrix")
    Threads.@threads for j in eachindex(peps)
        M[j, :] .= peptide_counts(kct, peps[j]; il_ambiguous=il_ambiguous, anchor=anchor)
        tick!(p)
    end
    return M
end
