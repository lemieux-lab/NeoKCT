using JSON

# Parses Jellyfish's own binary dump format (NOT this project's `JelloFish.jl` counter).
# See ARCHITECTURE.md for the on-disk layout and the packed-record bit math. Before trusting
# an unfamiliar Jellyfish build/version, spot-check with `verify_jellyfish_dump` below.

@inline function _read_le_u64(buf::AbstractVector{UInt8}, o::Int)::UInt64
    v = UInt64(0)
    @inbounds for i in 0:7
        v |= UInt64(buf[o+i]) << (8i)
    end
    return v
end

@inline function _read_le_uint(buf::AbstractVector{UInt8}, o::Int, n::Int)::UInt64
    v = UInt64(0)
    @inbounds for i in 0:n-1
        v |= UInt64(buf[o+i]) << (8i)
    end
    return v
end

# Zeroes the high padding bits so the raw word matches a real Kmer's right-aligned `.data[1]`
# (needed since it's used as a Dict key and later sorted/compared as a raw UInt64).
@inline function _clean_kmer_bits(bits::UInt64, K::Int, symbol_size::Int)::UInt64
    used = K * symbol_size
    mask = used >= 64 ? typemax(UInt64) : (UInt64(1) << used) - UInt64(1)
    return bits & mask
end

function _cmdline_int(header::AbstractDict, flags::String...)
    cmd = get(header, "cmdline", nothing)
    (isnothing(cmd) || isempty(cmd)) && return nothing
    for flag in flags
        i = findfirst(==(flag), cmd)
        isnothing(i) && continue
        i == length(cmd) && error("Jellyfish header cmdline flag $flag is missing its value")
        v = tryparse(Int, cmd[i+1])
        isnothing(v) && error("Jellyfish header cmdline flag $flag has a non-integer value: $(cmd[i+1])")
        return v
    end
    return nothing
end

function _is_canonical_dump(header::AbstractDict)
    c = get(header, "canonical", nothing)
    !isnothing(c) && return c === true || c == "true"
    cmd = get(header, "cmdline", String[])
    return "-C" in cmd || "--canonical" in cmd
end

"""
    read_jf_header(io::IO) -> (header::Dict, K::Int, count_bytes::Int)

Reads the leading `<offset>{...JSON...}` header of a Jellyfish binary dump from
`io`, leaving `io` positioned at the start of the first record. `K` is the
k-mer length (`-m`/`--mer-len`) and `count_bytes` the per-record counter width
(`--out-counter-len`, defaulting to Jellyfish's own default of 4) recovered
from the embedded `cmdline`.
"""
function read_jf_header(io::IO)
    offset_str = readuntil(io, "{")
    offset = tryparse(Int, offset_str)
    isnothing(offset) && error("Not a Jellyfish binary dump: expected an integer byte-offset before '{' " *
                                "at the start of the file")
    header_start = position(io) - 1  # position of the '{' consumed by readuntil
    seek(io, header_start)

    # Read as a String rather than JSON.parse(::IO): the latter errors on trailing
    # non-whitespace, which the binary record data right after the header would trip.
    header_bytes = Vector{UInt8}(undef, offset)
    read!(io, header_bytes)
    header_str = rstrip(c -> c == '\0' || isspace(c), String(header_bytes))
    header = JSON.parse(header_str)

    K = _cmdline_int(header, "-m", "--mer-len")
    isnothing(K) && error("Jellyfish header cmdline is missing -m/--mer-len (k-mer length)")
    count_bytes = something(_cmdline_int(header, "--out-counter-len"), 4)
    count_bytes in (1, 2, 4, 8) ||
        error("Unsupported Jellyfish counter length: $count_bytes byte(s) (expected 1, 2, 4 or 8)")
    return header, K, count_bytes
end

"""
    parse_jellyfish_dna_dump(fn; allow_canonical=false, min_count=0, record_chunk=1_000_000) -> (K, Dict{UInt64,UInt32})

Parses a Jellyfish binary dump into its k-mer length `K` and a `Dict` mapping
each DNA `K`-mer's raw `DNAAlphabet{2}` bit-encoding (right-aligned in the
`UInt64`, matching `Kmer{DNAAlphabet{2}, K, 1}.data[1]`) to its count
(saturating at `typemax(UInt32)`). Only `K <= 32` is supported, since a
`Kmer{DNAAlphabet{2}, K, 1}` must fit in a single `UInt64`.

`min_count` drops records with a count below the threshold (matches the old
`big_only` filter, generalised to an arbitrary cutoff). `record_chunk` controls
how many records are read from disk per bulk read (memory/syscall trade-off).

Errors if the dump was built with canonical counting (`-C`/`--canonical`)
unless `allow_canonical=true`. See `jellyfish_dump_hash` for why that matters
when the goal is single-frame translation to amino-acid k-mers.
"""
function parse_jellyfish_dna_dump(fn::String; allow_canonical::Bool=false, min_count::Integer=0,
                                   record_chunk::Int=1_000_000)
    isfile(fn) || error("No such file: $fn")
    record_chunk > 0 || throw(ArgumentError("record_chunk must be positive, got $record_chunk"))
    min_count >= 0 || throw(ArgumentError("min_count must be non-negative, got $min_count"))

    return open(fn, "r") do f
        header, K, count_bytes = read_jf_header(f)
        _is_canonical_dump(header) && !allow_canonical &&
            error("Jellyfish dump $fn was built with canonical counting (-C/--canonical): a canonical " *
                  "k-mer's stored bit pattern may be either strand, so translating it in a single fixed " *
                  "frame can silently yield the wrong peptide for roughly half the k-mers. Pass " *
                  "allow_canonical=true to proceed anyway.")

        symbol_size = Int(bits_per_symbol(DNAAlphabet{2}()))
        K * symbol_size <= 64 ||
            error("k-mer length K=$K ($(K*symbol_size) bits) exceeds 64 bits; k-mers longer than 32 " *
                  "bases are not supported by this parser")

        count_shift = 64 - count_bytes * 8
        count_mask = count_bytes == 8 ? typemax(UInt64) : (UInt64(1) << (count_bytes * 8)) - UInt64(1)
        packed = count_shift - K * symbol_size >= 0  # kmer + count share a single UInt64 word on disk
        record_size = packed ? 8 : 8 + count_bytes

        data_start = position(f)
        est = max(0, div(filesize(fn) - data_start, record_size))
        dna_counts = Dict{UInt64,UInt32}()
        sizehint!(dna_counts, est)

        min_count_u64 = UInt64(min_count)
        buf = Vector{UInt8}(undef, record_chunk * record_size)
        progress = ProgressUnknown(desc="Parsing Jellyfish dump $(basename(fn))...")
        while !eof(f)
            nb = readbytes!(f, buf, length(buf))
            nb == 0 && break
            n_full = nb ÷ record_size
            leftover = nb - n_full * record_size
            if leftover != 0
                eof(f) || error("Jellyfish dump $fn: misaligned record stream (got $nb bytes, not a " *
                                 "multiple of the $record_size-byte record size)")
                @warn "Jellyfish dump $fn: ignoring $leftover trailing byte(s) at EOF (truncated final record)"
            end
            @inbounds for r in 0:n_full-1
                o = r * record_size + 1
                raw = _read_le_u64(buf, o)
                count_raw = packed ? ((raw >> count_shift) & count_mask) : _read_le_uint(buf, o + 8, count_bytes)
                count_raw < min_count_u64 && continue
                kmer_bits = _clean_kmer_bits(raw, K, symbol_size)
                dna_counts[kmer_bits] = UInt32(min(count_raw, UInt64(typemax(UInt32))))
            end
            next!(progress; showvalues=[("unique k-mers seen", length(dna_counts))])
        end
        finish!(progress)
        return K, dna_counts
    end
end

"""
    jellyfish_dump_hash(fn; allow_canonical=false, min_count=0, record_chunk=1_000_000) -> Dict{UInt64,UInt32}

Parses a Jellyfish binary dump and translates every DNA k-mer to its
amino-acid k-mer using the same single-frame `translate` as
`jello_superthreaded_hash`/`count_kmers`, summing counts (saturating at
`typemax(UInt32)`) whenever more than one DNA k-mer (synonymous codons) maps
to the same AA k-mer. K-mers whose translation hits an in-frame stop codon are
dropped, same as the live-sequencing path.

The result has exactly the shape `KCT{K÷3, AAAlphabet}(...)` and
`push!(kct, ...)` expect, so a Jellyfish dump can be plugged directly into the
existing KCT-building pipeline in place of a `jello_superthreaded_hash` sample.
Useful when only pre-computed counts (e.g. GTEx) are available and the raw
reads aren't.

See `parse_jellyfish_dna_dump` for the `allow_canonical`/`min_count`/
`record_chunk` keywords, and its docstring for why canonical counting is
rejected by default.
"""
function jellyfish_dump_hash(fn::String; allow_canonical::Bool=false, min_count::Integer=0,
                              record_chunk::Int=1_000_000)
    K, dna_counts = parse_jellyfish_dna_dump(fn; allow_canonical, min_count, record_chunk)

    aa_counts = Dict{UInt64,UInt32}()
    sizehint!(aa_counts, length(dna_counts))
    n_stop = 0
    @showprogress desc="Translating Jellyfish DNA $K-mers to AA $(K÷3)-mers..." for (kbits, cnt) in dna_counts
        kmer = Kmer{DNAAlphabet{2}, K, 1}(Kmers.unsafe, (kbits,))
        aa_kmer = translate(kmer)
        if isnothing(aa_kmer)
            n_stop += 1
            continue
        end
        aa_bits = aa_kmer.data[1]
        prev = get(aa_counts, aa_bits, UInt32(0))
        aa_counts[aa_bits] = UInt32(min(UInt64(prev) + UInt64(cnt), UInt64(typemax(UInt32))))
    end

    printstyled("Jellyfish dump $fn: $(length(dna_counts)) DNA $K-mers -> $(length(aa_counts)) unique " *
                "AA $(K÷3)-mers ($n_stop dropped for in-frame stop codons)\n", color=:green)
    return aa_counts
end

"""
    verify_jellyfish_dump(fn; n::Int=1000) -> Bool

Best-effort sanity check: parses the first `n` records with
`parse_jellyfish_dna_dump` and compares them against `jellyfish dump -c fn`
(the k-mer/count columns, text format) run via the `jellyfish` CLI, if it's on
`PATH`. Returns `true` on a full match, `false` on any mismatch (printing the
first disagreement), and `missing` if the `jellyfish` binary isn't available
to check against.

Useful before trusting a full production parse against an unfamiliar
Jellyfish build/version, since the binary record layout this parser assumes
was reverse-engineered rather than read from a format spec.
"""
function verify_jellyfish_dump(fn::String; n::Int=1000)
    isnothing(Sys.which("jellyfish")) && return missing

    K, parsed = parse_jellyfish_dna_dump(fn; allow_canonical=true, record_chunk=max(n, 1))
    expected = Dict{UInt64,UInt32}()
    open(`jellyfish dump -c -L 0 $fn`) do io
        for (i, line) in enumerate(eachline(io))
            i > n && break
            fields = split(line)
            length(fields) == 2 || error("Unexpected `jellyfish dump -c` line: $line")
            seq, count = fields
            length(seq) == K || error("Unexpected k-mer length in `jellyfish dump -c` output: $seq")
            bits = Kmer{DNAAlphabet{2}, K}(LongDNA{2}(seq)).data[1]
            expected[bits] = UInt32(parse(Int, count))
        end
    end

    for (bits, count) in expected
        got = get(parsed, bits, nothing)
        if got != count
            printstyled("Mismatch at k-mer bits $bits: expected count $count, parsed $got\n", color=:red)
            return false
        end
    end
    return true
end
