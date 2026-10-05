## Sample metadata (library sizes for rphm normalization, manuscript Section 2.7) ##

# Per-sample tissue label and normalization denominator, in table slot order. `lib` is
# BamQuery's own denominator (total primary alignments, from A1's flagstat pass), so rphm
# values are directly comparable to BamQuery's own 8.55 rphm threshold.
struct SampleMeta
    sample_id::Vector{String}
    tissue::Vector{String}
    lib::Vector{Float64}
end

"""
    read(io::IO, ::Type{SampleMeta}; lib_col="primary_reads") -> SampleMeta

Read a tab-separated sample manifest (`#`-comment lines and blank lines skipped) with at
least `slot`, `sample_id`, `tissue` and `lib_col` columns, sorted into table slot order.
A blank `lib_col` value (a sample whose library size is unknown, see A1) reads as `NaN`,
not zero, so it can't silently pass or fail a threshold: callers should drop `NaN` rows
before aggregating rather than propagate them into an rphm value.
"""
function Base.read(io::IO, ::Type{SampleMeta}; lib_col::String="primary_reads")
    # split on the un-stripped line: stripping first would eat a trailing empty field (a
    # sample with no known library size, see A1) along with the newline, desyncing columns.
    rows = [split(chomp(l), '\t') for l in eachline(io) if !startswith(l, '#') && !isempty(chomp(l))]
    hdr = String.(rows[1])
    col(n) = findfirst(==(n), hdr)
    body = rows[2:end]
    sort!(body; by=r -> parse(Int, r[col("slot")]))
    lib = [isempty(v) ? NaN : parse(Float64, v) for v in getindex.(body, col(lib_col))]
    return SampleMeta(String.(getindex.(body, col("sample_id"))),
                      String.(getindex.(body, col("tissue"))), lib)
end

Base.length(meta::SampleMeta) = length(meta.sample_id)

"""
    rphm(M::AbstractMatrix{<:Integer}, meta::SampleMeta) -> Matrix{Float64}

Reads per hundred million, BamQuery's normalized unit: each column (sample) of `M`
(peptides x samples, raw counts) divided by that sample's library size and scaled by 1e8.
A `NaN` library size (see `read`) propagates to `NaN` in that column, not a false zero.
"""
rphm(M::AbstractMatrix{<:Integer}, meta::SampleMeta) = M ./ reshape(meta.lib, 1, :) .* 1e8
