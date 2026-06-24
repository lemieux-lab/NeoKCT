if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(".")
end

using GZip
using EzXML

abstract type Fq end

struct FastqGZ <: Fq
    iterable::GZipStream
end

struct Fastq <: Fq
    iterable::IOStream
end

struct Mzid
    iterable::Vector{EzXML.Node}
end

function Base.iterate(seqfile::F, state=nothing) where {F<:Fq}
    eof(seqfile.iterable) && return nothing
    readline(seqfile.iterable)
    seq = strip(readline(seqfile.iterable))
    sep = readline(seqfile.iterable)
    startswith(sep, '+') || error("FASTQ parser desync: expected '+' separator, got $(repr(sep)); sequence was $(repr(seq))")
    readline(seqfile.iterable)
    return seq, state
end

function Base.iterate(seqfile::Mzid, state::Int=1)
    return state < length(seqfile.iterable) ? (seqfile.iterable[state].content, state+1) : nothing
end

stream(::Type{Val{:fastq}}, path::String) = Fastq(open(path, "r"))
stream(::Type{Val{:fq}}, path::String) = Fastq(open(path, "r"))
stream(::Type{Val{:mzid}}, path::String) = Mzid(findall("//*[local-name()='Seq']", readxml(path)))
stream(::Type{Val{:fa}}, path::String) = FastaFile(open(path, "r"))
stream(::Type{Val{:fasta}}, path::String) = FastaFile(open(path, "r"))
stream(::Type{Val{:gtf}}, path::String) = GtfFile(open(path, "r"))
stream(::Type{Val{:gff}}, path::String) = Gff3File(open(path, "r"))
stream(::Type{Val{:gff3}}, path::String) = Gff3File(open(path, "r"))

function stream(path::String)
    parts = split(path, ".")
    ext = parts[end]
    if ext == "gz" && length(parts) >= 2
        inner = parts[end-1]
        (inner == "fastq" || inner == "fq") && return FastqGZ(gzopen(path, "r"))
        (inner == "fa" || inner == "fasta") && return FastaFileGZ(gzopen(path, "r"))
        inner == "gtf" && return GtfFileGZ(gzopen(path, "r"))
        (inner == "gff" || inner == "gff3") && return Gff3FileGZ(gzopen(path, "r"))
        return FastqGZ(gzopen(path, "r"))  # legacy fallback
    end
    return stream(Val{Symbol(ext)}, path)
end

# ---- FASTA ----

struct FastaRecord
    header::String
    sequence::String
end

struct FastaFile
    iterable::IOStream
end

struct FastaFileGZ
    iterable::GZipStream
end

function _iterate_fasta(io, state::Union{Nothing, String})
    header = if isnothing(state)
        h = nothing
        while !eof(io)
            line = strip(readline(io))
            isempty(line) && continue
            if startswith(line, '>')
                h = String(line[2:end]); break
            end
        end
        isnothing(h) && return nothing
        h
    else
        state
    end
    seq = IOBuffer()
    next_header = nothing
    while !eof(io)
        line = strip(readline(io))
        isempty(line) && continue
        if startswith(line, '>')
            next_header = String(line[2:end]); break
        end
        write(seq, line)
    end
    return FastaRecord(header, String(take!(seq))), next_header
end

Base.iterate(fa::FastaFile, state::Union{Nothing, String}=nothing) = _iterate_fasta(fa.iterable, state)
Base.iterate(fa::FastaFileGZ, state::Union{Nothing, String}=nothing) = _iterate_fasta(fa.iterable, state)
Base.IteratorSize(::Type{FastaFile}) = Base.SizeUnknown()
Base.IteratorSize(::Type{FastaFileGZ}) = Base.SizeUnknown()

# ---- GTF ----

struct GtfRecord
    seqname::String
    feature::String
    start::Int
    stop::Int
    strand::Char
    attrs::Dict{String, String}
end

struct GtfFile
    iterable::IOStream
end

struct GtfFileGZ
    iterable::GZipStream
end

function _parse_gtf_attrs(s::AbstractString)
    attrs = Dict{String, String}()
    for m in eachmatch(r"(\w+)\s+\"([^\"]+)\"", s)
        attrs[m.captures[1]] = m.captures[2]
    end
    return attrs
end

function _iterate_gtf(io, state)
    while !eof(io)
        line = readline(io)
        startswith(line, '#') && continue
        isempty(strip(line)) && continue
        f = split(line, '\t')
        length(f) < 9 && continue
        return GtfRecord(f[1], f[3], parse(Int, f[4]), parse(Int, f[5]),
                         f[7][1], _parse_gtf_attrs(f[9])), state
    end
    return nothing
end

Base.iterate(gtf::GtfFile, state=nothing) = _iterate_gtf(gtf.iterable, state)
Base.iterate(gtf::GtfFileGZ, state=nothing) = _iterate_gtf(gtf.iterable, state)
Base.IteratorSize(::Type{GtfFile}) = Base.SizeUnknown()
Base.IteratorSize(::Type{GtfFileGZ}) = Base.SizeUnknown()

# ---- GFF3 ----

struct Gff3Record
    seqname::String
    feature::String
    start::Int
    stop::Int
    strand::Char
    attrs::Dict{String, String}
end

struct Gff3File
    iterable::IOStream
end

struct Gff3FileGZ
    iterable::GZipStream
end

function _parse_gff3_attrs(s::AbstractString)
    attrs = Dict{String, String}()
    for kv in split(s, ';')
        kv = strip(kv)
        isempty(kv) && continue
        eq = findfirst('=', kv)
        isnothing(eq) && continue
        attrs[String(kv[1:eq-1])] = String(kv[eq+1:end])
    end
    return attrs
end

function _iterate_gff3(io, state)
    while !eof(io)
        line = readline(io)
        startswith(line, '#') && continue
        strip(line) == "##FASTA" && return nothing  # embedded FASTA section ends records
        isempty(strip(line)) && continue
        f = split(line, '\t')
        length(f) < 9 && continue
        return Gff3Record(f[1], f[3], parse(Int, f[4]), parse(Int, f[5]),
                          f[7][1], _parse_gff3_attrs(f[9])), state
    end
    return nothing
end

Base.iterate(gff::Gff3File, state=nothing) = _iterate_gff3(gff.iterable, state)
Base.iterate(gff::Gff3FileGZ, state=nothing) = _iterate_gff3(gff.iterable, state)
Base.IteratorSize(::Type{Gff3File}) = Base.SizeUnknown()
Base.IteratorSize(::Type{Gff3FileGZ}) = Base.SizeUnknown()
