## Sparse Counts Layer (.kct V4.0) ##

"""
    SparseCountsLayer <: AbstractLayer

Row-deduplicated, sparse store of per-k-mer count vectors (`.kct` V4.0).

Each k-mer maps, via `row_id`, to one row in a shared pool. A row is the list of
`(sample, count)` pairs for the samples where that k-mer is present, in ascending sample
order. Absent samples are implicit and never stored. K-mers whose count vectors are
identical share a single pooled row.

Unlike `CountsLayer` (V3.0, `flat_cids` + `n_cids` + a per-k-mer `PackedArray` word chain
with materialised absent-sample zeros), nothing here grows with the sample count except the
number of pairs inside a row. `row_id` is one `UInt32` per k-mer, and the pool grows only
with the number of distinct count vectors. Built only by `build_kct_streaming`. A
`SparseCountsLayer` KCT is read-only (`push!` / `collapse!` / `repack` / `sort!` all error).

# Fields
- `row_id`: 1-based pool index per k-mer, `length(row_id)` is the k-mer count
- `row_offsets`: CSR offsets into `row_samples` / `row_counts`, `length` is `n_rows + 1`,
  1-based, `row_offsets[1] == 1`
- `row_samples`: flat pool of sample indices, ascending within each row
- `row_counts`: flat pool of counts, parallel to `row_samples`
- `n_samples`: total number of samples represented (used to size the zero-padded read)
"""
mutable struct SparseCountsLayer <: AbstractLayer
    row_id::Vector{UInt32}
    row_offsets::Vector{UInt64}
    row_samples::Vector{UInt32}
    row_counts::Vector{UInt32}
    n_samples::Base.RefValue{Int64}
end

SparseCountsLayer(row_id::Vector{UInt32}, row_offsets::Vector{UInt64},
                  row_samples::Vector{UInt32}, row_counts::Vector{UInt32}, n_samples::Integer) =
    SparseCountsLayer(row_id, row_offsets, row_samples, row_counts, Ref(Int64(n_samples)))

Base.length(scl::SparseCountsLayer) = length(scl.row_id)

# Number of distinct pooled rows.
_n_rows(scl::SparseCountsLayer) = length(scl.row_offsets) - 1

# Full length-`n_samples` count vector for k-mer `i`: scatter its pooled row's
# `(sample, count)` pairs into a zero vector. Same contract as the CountsLayer method in
# KCTLayers.jl, but always a `Vector{UInt32}`. V3.0 can promote to the packed word type via
# `vcat` when it materialises trailing zeros, while V4.0's plain `UInt32` matches the test
# literals and `_expected_cv`.
function assemble_count_vector(scl::SparseCountsLayer, i::Integer)
    r = scl.row_id[i]
    lo = scl.row_offsets[r]
    hi = scl.row_offsets[r + 1] - 1
    v = zeros(UInt32, scl.n_samples.x)
    @inbounds for k in lo:hi
        v[scl.row_samples[k]] = scl.row_counts[k]
    end
    return v
end

Base.getindex(scl::SparseCountsLayer, i::Integer) = assemble_count_vector(scl, i)
