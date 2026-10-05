# Builds a GenomicIndex from an Ensembl transcript FASTA + GTF or GFF3 annotation.
#
# Bitmask layout: bit b-1 is set when a k-mer comes from biotype_names[b].
# biotype_names[1] is always "intergenic" (INTERGENIC_MASK = 1 << 0 = 1).
# k-mers appearing in multiple biotypes accumulate bits via OR.

# GTF: transcript lines carry a `transcript_biotype` attribute.
# GFF3: transcript/mRNA lines carry a `biotype` attribute, with an ID like "transcript:ENST...".
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
    build_genomic_index(fasta_paths, annotation_path, K; proteome_path=nothing, ...)

Build a `KCT{K÷3, AAAlphabet, Nothing, BiotypLayer}` from one or more Ensembl transcript
FASTAs (a single path, or a vector -- e.g. `cdna.all.fa` and `ncrna.fa` together, so
lncRNA/other ncRNA transcripts stop defaulting to intergenic) and a GTF or GFF3.
`K` is the nucleotide k-mer length, same convention as `build_kct`.

Biotype membership comes from `transcript_biotype` in the FASTA headers (Ensembl cdna.all.fa),
falling back to the annotation file. A k-mer that appears in more than one biotype accumulates
its bits via OR. A k-mer with no annotation defaults to `INTERGENIC_MASK`.

`proteome_path`, if given (Ensembl `pep.all.fa`), adds an `"in_frame_cds"` bit: every
`K÷3`-residue window of every annotated protein is by definition an in-frame canonical
k-mer, giving the canonical-vs-out-of-frame split without parsing CDS coordinates.
"""
function build_genomic_index(fasta_paths::Union{String, Vector{String}}, annotation_path::String, K::Int;
                              proteome_path::Union{String, Nothing}=nothing,
                              checkpoint_size::Type{<:Unsigned}=UInt64,
                              delta_size::Type{<:Unsigned}=UInt32)
    sorted_kmers, bitmasks, biotype_names = _build_genomic_index_data(fasta_paths, annotation_path, K;
                                                                      proteome_path=proteome_path)
    return KCT{K÷3, AAAlphabet}(sorted_kmers, bitmasks, biotype_names;
                                checkpoint_size=checkpoint_size,
                                delta_size=delta_size)
end

# Sorted, unique AA K-mer codes of every window of the reference proteome -- every one is by
# definition in-frame canonical CDS, no need to parse CDS coordinates off the transcript.
function proteome_kmers(pep_fasta::String, K::Int)
    ks = UInt64[]
    prog = ProgressUnknown(desc="Scanning reference proteome for in_frame_cds k-mers...")
    for rec in stream(pep_fasta)
        rec isa FastaRecord || continue
        seq = rec.sequence
        for s in 1:(length(seq) - K + 1)
            w = seq[s:s + K - 1]
            occursin(r"[^ACDEFGHIKLMNPQRSTVWY]", w) && continue
            push!(ks, _aa_code(w))
        end
        next!(prog)
    end
    finish!(prog)
    return unique!(sort!(ks))
end

function _build_genomic_index_data(fasta_paths::Union{String, Vector{String}}, annotation_path::String, K::Int;
                                    proteome_path::Union{String, Nothing}=nothing)
    paths = fasta_paths isa String ? [fasta_paths] : fasta_paths
    tid_to_biotype = _parse_transcript_biotypes(annotation_path)

    # Collect all encountered biotype names, with intergenic first. Pre-scan to build a
    # stable ordered list before assigning bit positions. "in_frame_cds" is appended last
    # (not sorted in with the rest) so its bit position doesn't depend on which transcript
    # biotypes happen to be present in a given annotation release.
    seen_biotypes = Set{String}()
    for (_, bt) in tid_to_biotype; push!(seen_biotypes, bt); end
    biotype_names = vcat(["intergenic"], sort!(collect(seen_biotypes)))
    !isnothing(proteome_path) && push!(biotype_names, "in_frame_cds")
    biotype_index = Dict(name => i for (i, name) in enumerate(biotype_names))  # 1-based

    kmer_masks = Dict{UInt64, UInt64}()

    for fasta_path in paths
        progress = ProgressUnknown(desc="Building GenomicIndex k-mer table ($(basename(fasta_path)))...")
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
    end

    if !isnothing(proteome_path)
        in_frame_bit = UInt64(1) << (biotype_index["in_frame_cds"] - 1)
        for k_bits in proteome_kmers(proteome_path, K ÷ 3)
            kmer_masks[k_bits] = get(kmer_masks, k_bits, UInt64(0)) | in_frame_bit
        end
    end

    sorted_kmers = sort!(collect(keys(kmer_masks)))
    bitmasks = [kmer_masks[k] for k in sorted_kmers]
    printstyled("GenomicIndex data: $(length(sorted_kmers)) unique k-mers, $(length(biotype_names)) biotypes\n",
                color=:green)
    return sorted_kmers, bitmasks, biotype_names
end
