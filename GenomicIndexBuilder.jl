# Builds a GenomicIndex from an Ensembl transcript FASTA + GTF or GFF3 annotation.
#
# Bitmask layout: bit b-1 is set when a k-mer comes from biotype_names[b].
# biotype_names[1] is always "intergenic" (INTERGENIC_MASK = 1 << 0 = 1).
# k-mers appearing in multiple biotypes accumulate bits via OR.

# GTF: transcript lines carry a `transcript_biotype` attribute.
# GFF3: transcript/mRNA lines carry a `biotype` attribute; ID looks like "transcript:ENST...".
function _parse_transcript_biotypes_gtf(path::String)
    tid_to_biotype = Dict{String, String}()
    for rec in stream(path)
        rec isa GtfRecord || continue
        rec.feature == "transcript" || continue
        tid = get(rec.attrs, "transcript_id", nothing)
        bt = get(rec.attrs, "transcript_biotype", nothing)
        (isnothing(tid) || isnothing(bt)) && continue
        tid_to_biotype[tid] = bt
    end
    return tid_to_biotype
end

function _parse_transcript_biotypes_gff3(path::String)
    tid_to_biotype = Dict{String, String}()
    for rec in stream(path)
        rec isa Gff3Record || continue
        (rec.feature == "transcript" || rec.feature == "mRNA") || continue
        raw_id = get(rec.attrs, "ID", nothing)
        bt = get(rec.attrs, "biotype", nothing)
        (isnothing(raw_id) || isnothing(bt)) && continue
        # Ensembl GFF3 IDs are prefixed: "transcript:ENST..."
        tid = startswith(raw_id, "transcript:") ? raw_id[12:end] : raw_id
        tid_to_biotype[tid] = bt
    end
    return tid_to_biotype
end

function _parse_transcript_biotypes(annotation_path::String)
    ext = split(annotation_path, ".")[end]
    inner = ext == "gz" ? split(annotation_path, ".")[end-1] : ext
    inner == "gtf" && return _parse_transcript_biotypes_gtf(annotation_path)
    (inner == "gff" || inner == "gff3") && return _parse_transcript_biotypes_gff3(annotation_path)
    error("Unsupported annotation format: $annotation_path")
end

# Extract transcript ID from Ensembl FASTA header:
# ">ENST00000641515.2 cdna chromosome:... transcript_biotype:protein_coding"
# Falls back to GTF biotype map if no transcript_biotype in header.
function _header_biotype(header::String)
    m = match(r"transcript_biotype:(\S+)", header)
    isnothing(m) ? nothing : m.captures[1]
end

function _header_tid(header::String)
    tid = split(header, ' ')[1]
    # Strip version suffix (ENST00000641515.2 → ENST00000641515)
    dot = findfirst('.', tid)
    isnothing(dot) ? tid : tid[1:dot-1]
end

"""
    build_genomic_index(fasta_path, annotation_path, K; ...)

Build a `GenomicIndex{K÷3, AAAlphabet}` from an Ensembl transcript FASTA and a GTF or GFF3.
`K` is the **DNA** k-mer length (same convention as `build_kct`).

Biotype membership is determined from `transcript_biotype` in the FASTA headers (Ensembl cdna.all.fa)
or, as a fallback, from the annotation file. k-mers appearing in multiple biotypes accumulate bits
via OR. k-mers with no annotation default to `INTERGENIC_MASK`.
"""
function build_genomic_index(fasta_path::String, annotation_path::String, K::Int;
                              checkpoint_size::Type{<:Unsigned}=UInt64,
                              delta_size::Type{<:Unsigned}=UInt32)
    sorted_kmers, bitmasks, biotype_names = _build_genomic_index_data(fasta_path, annotation_path, K)
    return GenomicIndex{K÷3, AAAlphabet}(sorted_kmers, bitmasks, biotype_names;
                                          checkpoint_size=checkpoint_size,
                                          delta_size=delta_size)
end

function _build_genomic_index_data(fasta_path::String, annotation_path::String, K::Int)
    tid_to_biotype = _parse_transcript_biotypes(annotation_path)

    # Collect all encountered biotype names; intergenic must be first.
    # We pre-scan to build a stable ordered list before assigning bit positions.
    seen_biotypes = Set{String}()
    for (_, bt) in tid_to_biotype; push!(seen_biotypes, bt); end
    biotype_names = vcat(["intergenic"], sort!(collect(seen_biotypes)))
    biotype_index = Dict(name => i for (i, name) in enumerate(biotype_names))  # 1-based

    kmer_masks = Dict{UInt64, UInt64}()

    progress = ProgressUnknown(desc="Building GenomicIndex k-mer table...")
    for rec in stream(fasta_path)
        rec isa FastaRecord || continue
        isempty(rec.sequence) && continue
        'N' in rec.sequence && continue

        biotype = _header_biotype(rec.header)
        if isnothing(biotype)
            tid = _header_tid(rec.header)
            biotype = get(tid_to_biotype, tid, "intergenic")
        end
        bit_pos = get(biotype_index, biotype, 1)  # fallback to intergenic if unknown
        bit = UInt64(1) << (bit_pos - 1)

        seq = LongSequence{DNAAlphabet{2}}(rec.sequence)
        length(seq) < K && continue
        for kmer in k_merize(seq, K=K)
            aa_kmer = translate(kmer)
            isnothing(aa_kmer) && continue
            k_bits = aa_kmer.data[1]
            kmer_masks[k_bits] = get(kmer_masks, k_bits, UInt64(0)) | bit
        end
        next!(progress)
    end
    finish!(progress)

    sorted_kmers = sort!(collect(keys(kmer_masks)))
    bitmasks = [kmer_masks[k] for k in sorted_kmers]
    printstyled("GenomicIndex data: $(length(sorted_kmers)) unique k-mers, $(length(biotype_names)) biotypes\n",
                color=:green)
    return sorted_kmers, bitmasks, biotype_names
end
