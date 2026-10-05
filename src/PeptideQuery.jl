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

"""
    peptide_counts(kct, pep; il_ambiguous=false) -> Vector{UInt32}

Per-sample abundance of `pep` against a `CountsLayer` table. `length(pep) <= K` sums the
counts over the peptide's whole prefix range: the number of read windows that start with a
coding sequence of `pep` and continue stop-free for `3(K - length(pep))` more nucleotides,
mirroring BamQuery's requirement that a read span the whole MCS. `length(pep) == K` is the
k-mer's own row. `length(pep) > K` takes the elementwise minimum over the
`length(pep) - K + 1` overlapping K-mer tiles, an upper bound on the number of read windows
carrying the whole peptide. `il_ambiguous=true` sums over every I/L variant of `pep`, for a
peptide whose I/L assignment came from a database lookup rather than the spectrum (leucine
and isoleucine are isobaric, indistinguishable by mass alone).

Two known blind spots, not fixed here (PAPER_TODO.md task A5): a peptide sitting right at
the C-terminus of its ORF has no stop-free extension and returns 0 for `length(pep) < K`
even if genuinely present; matching is exact, so a read with a sequencing error inside the
window is not counted, same as BamQuery.
"""
function peptide_counts(kct::KCT{K, AAAlphabet, CountsLayer}, pep::AbstractString;
                        il_ambiguous::Bool=false) where {K}
    il_ambiguous && return sum(peptide_counts(kct, v) for v in _il_variants(pep))
    L = length(pep)
    L <= K && return _sum_rows(kct.counts, searchsorted(kct.kmer, pep))
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
    peptide_matrix(kct, peps; il_ambiguous=false) -> Matrix{UInt32}

`peptide_counts` for every peptide in `peps`, one row per peptide, one column per sample.
Peptides are independent, so the loop runs across threads. Logs progress through `_Prog`
(Progress.jl), since a real peptide set against the full GTEx table is exactly the kind of
run PAPER_TODO.md's conventions want instrumented.
"""
function peptide_matrix(kct::KCT{K, AAAlphabet, CountsLayer}, peps::Vector{String};
                        il_ambiguous::Bool=false) where {K}
    M = zeros(UInt32, length(peps), Int(kct.counts.n_samples.x))
    p = _Prog(length(peps), "peptide_matrix")
    Threads.@threads for j in eachindex(peps)
        M[j, :] .= peptide_counts(kct, peps[j]; il_ambiguous=il_ambiguous)
        tick!(p)
    end
    return M
end
