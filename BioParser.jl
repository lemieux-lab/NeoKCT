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
    readline(seqfile.iterable)
    readline(seqfile.iterable)
    return seq, state
end

function Base.iterate(seqfile::Mzid, state::Int=1)
    return state < length(seqfile.iterable) ? (seqfile.iterable[state].content, state+1) : nothing
end

stream(::Type{Val{:fastq}}, path::String) = Fastq(open(path, "r"))
stream(::Type{Val{:gz}}, path::String) = FastqGZ(gzopen(path, "r"))  # TODO: Assumes fastq.gz
stream(::Type{Val{:mzid}}, path::String) = Mzid(findall("//*[local-name()='Seq']", readxml(path)))
stream(path::String) = stream(Val{Symbol(split(path, ".")[end])}, path)
