using Test
using Kmers, BioSequences, BioSymbols, NArrays, ProgressMeter, JSON

include("../src/KCTLayers.jl")

# Builds a minimal synthetic Jellyfish binary/sorted dump: <offset>{json header}<records>,
# matching what `jellyfish count`/`jellyfish dump` (without -c/--fasta) write. `records` is a
# list of (DNA::String, count::Integer) pairs. Mirrors parse_jellyfish_dna_dump's own packed-vs-
# separate branch condition so callers don't need to reason about which layout a given K/
# count_bytes combination selects.
function _write_jf_dump(path::String, K::Int, count_bytes::Int, records; canonical::Bool=false)
    open(path, "w") do io
        cmdline = ["jellyfish", "count", "-m", string(K), "-o", "out.jf",
                   "--out-counter-len", string(count_bytes), "-s", "100M"]
        canonical && push!(cmdline, "-C")
        json_str = JSON.json(Dict("cmdline" => cmdline, "format" => "binary/sorted"))
        write(io, string(length(json_str)))
        write(io, json_str)

        symbol_size = Int(bits_per_symbol(DNAAlphabet{2}()))
        count_shift = 64 - count_bytes * 8
        packed = count_shift - K * symbol_size >= 0
        for (seq, count) in records
            kbits = Kmer{DNAAlphabet{2}, K}(LongDNA{2}(seq)).data[1]
            if packed
                write(io, kbits | (UInt64(count) << count_shift))
            else
                write(io, kbits)
                write(io, (count_bytes == 1 ? UInt8 : count_bytes == 2 ? UInt16 :
                           count_bytes == 4 ? UInt32 : UInt64)(count))
            end
        end
    end
end

# K=4 amino acid k-mers (20-bit encoding) keep the prefix index at 1 entry,
# which lets us test the full KCT stack without large memory allocations.
# The bit patterns below are arbitrary sorted UInt64 values < 2^20.

@testset "NeoKCT" begin

  @testset "1. _kmer_offsets invariants" begin
    n_cids = UInt16[1, 2, 0, 3, 1]
    offsets = _kmer_offsets(n_cids)

    @test offsets[1] == 1
    @test offsets[end] == 1 + sum(n_cids)
    for i in eachindex(n_cids)
      @test offsets[i+1] - offsets[i] == n_cids[i]
    end

    # empty input
    @test _kmer_offsets(UInt16[])[1] == 1
    @test length(_kmer_offsets(UInt16[])) == 1
  end

  @testset "2. BiotypLayer bitmask logic" begin
    bnames = ["intergenic", "protein_coding", "lncRNA"]
    bl = BiotypLayer(3, bnames)

    @test all(==(UInt16(1)), bl.ids)
    @test biotype_names_for(bl, 1) == ["intergenic"]
    @test has_biotype(bl, 1, "intergenic")
    @test !has_biotype(bl, 1, "protein_coding")
    @test_throws ArgumentError has_biotype(bl, 1, "unknown_biotype")

    # Intern a mask that combines protein_coding (bit 1) and lncRNA (bit 2).
    # biotype_names_for maps bit b-1 → biotype_names[b], so:
    #   protein_coding = b=2 → bit 1 → 0b010
    #   lncRNA         = b=3 → bit 2 → 0b100
    pool = copy(bl.pool)
    idx_map = Dict{UInt64, UInt16}(pool[i] => UInt16(i) for i in eachindex(pool))
    pc_lnc_mask = UInt64(0b110)
    id2 = _intern_mask!(pool, idx_map, pc_lnc_mask)

    bl2 = BiotypLayer([UInt16(1), id2, UInt16(1)], pool, bnames)
    @test sort(biotype_names_for(bl2, 2)) == sort(["protein_coding", "lncRNA"])
    @test has_biotype(bl2, 2, "protein_coding")
    @test has_biotype(bl2, 2, "lncRNA")
    @test !has_biotype(bl2, 2, "intergenic")

    # _intern_mask! is idempotent: same mask returns same id
    @test _intern_mask!(pool, idx_map, pc_lnc_mask) == id2
  end

  @testset "3. KCT single-sample build and lookup" begin
    raw = Dict{UInt64, UInt32}(
      UInt64(100) => UInt32(3),
      UInt64(200) => UInt32(7),
      UInt64(500) => UInt32(1),
      UInt64(1000) => UInt32(9),
    )
    kct = KCT{4, AAAlphabet}(raw)

    @test length(kct) == 4

    # k-mers must come out sorted
    @test issorted(collect(kct.kmer.seqs))

    # every inserted k-mer must be found at a valid position
    for (k_bits, count) in raw
      pos = findfirst(kct.kmer, k_bits)
      @test pos > 0
      cv = kct.counts[pos]
      @test cv[1] == count
    end

    # absent k-mer must return 0
    @test findfirst(kct.kmer, UInt64(999)) == 0
  end

  @testset "4. Two-sample push! and count reconstruction" begin
    s1 = Dict{UInt64, UInt32}(UInt64(100) => UInt32(3), UInt64(200) => UInt32(5))
    s2 = Dict{UInt64, UInt32}(UInt64(200) => UInt32(2), UInt64(500) => UInt32(7))

    kct = KCT{4, AAAlphabet}(s1)
    push!(kct, s2)

    @test length(kct) == 3
    @test kct.counts.samples[] == 2

    pos_100 = findfirst(kct.kmer, UInt64(100))
    pos_200 = findfirst(kct.kmer, UInt64(200))
    pos_500 = findfirst(kct.kmer, UInt64(500))
    @test pos_100 > 0 && pos_200 > 0 && pos_500 > 0

    cv_100 = kct.counts[pos_100]
    cv_200 = kct.counts[pos_200]
    cv_500 = kct.counts[pos_500]

    # 100: present s1=3, absent s2 → trailing zero
    @test cv_100[1] == 3
    @test cv_100[2] == 0

    # 200: present in both
    @test cv_200[1] == 5
    @test cv_200[2] == 2

    # 500: absent s1 → backfilled zero, present s2=7
    @test cv_500[1] == 0
    @test cv_500[2] == 7

    # sorted invariant must hold after merge
    @test issorted(collect(kct.kmer.seqs))
  end

  @testset "5. collapse! preserves counts and deduplicates" begin
    # Two k-mers with the same count get deduplicated to one PackedArray word.
    s1 = Dict{UInt64, UInt32}(UInt64(100) => UInt32(5), UInt64(200) => UInt32(5))
    kct = KCT{4, AAAlphabet}(s1)

    pos_100 = findfirst(kct.kmer, UInt64(100))
    pos_200 = findfirst(kct.kmer, UInt64(200))
    cv_before_100 = copy(kct.counts[pos_100])
    cv_before_200 = copy(kct.counts[pos_200])

    n_words_before = length(kct.counts.counts.words)
    kct2 = collapse!(kct)

    pos_100c = findfirst(kct2.kmer, UInt64(100))
    pos_200c = findfirst(kct2.kmer, UInt64(200))

    @test kct2.counts[pos_100c] == cv_before_100
    @test kct2.counts[pos_200c] == cv_before_200

    # identical count vectors must have been folded to fewer words
    @test length(kct2.counts.counts.words) < n_words_before
  end

  @testset "6. KCTLoader write/load round-trip" begin
    raw = Dict{UInt64, UInt32}(
      UInt64(100) => UInt32(3),
      UInt64(200) => UInt32(7),
      UInt64(500) => UInt32(1),
    )
    kct = KCT{4, AAAlphabet}(raw)

    path = tempname() * ".kct"
    try
      write_kct(kct, path)
      kct2 = load_kct(path)

      @test length(kct2) == length(kct)
      @test issorted(collect(kct2.kmer.seqs))

      for k_bits in keys(raw)
        pos1 = findfirst(kct.kmer, k_bits)
        pos2 = findfirst(kct2.kmer, k_bits)
        @test pos1 > 0 && pos2 > 0
        @test kct.counts[pos1] == kct2.counts[pos2]
      end
    finally
      isfile(path) && rm(path)
    end
  end

  @testset "7. JellyfishDump parsing" begin
    K = 6; count_bytes = 1  # K*2 + count_bytes*8 = 12+8 = 20 <= 64: packed single-word records

    @testset "packed-record round trip (kmer + count share one UInt64 word)" begin
      records = [("ATGCCC", 4), ("AAACCC", 200)]
      path = tempname()
      try
        _write_jf_dump(path, K, count_bytes, records)
        Kp, dna = parse_jellyfish_dna_dump(path)
        @test Kp == K
        for (seq, count) in records
          bits = Kmer{DNAAlphabet{2}, K}(LongDNA{2}(seq)).data[1]
          @test dna[bits] == UInt32(count)
        end
      finally
        isfile(path) && rm(path)
      end
    end

    @testset "separate-fields record round trip (K=31, count_bytes=4)" begin
      # K*2 + count_bytes*8 = 62+32 = 94 > 64: kmer occupies its own 8-byte word, count follows separately.
      records = [("A"^31, 3), ("ACGT"^7 * "AAA", 70_000)]
      path = tempname()
      try
        _write_jf_dump(path, 31, 4, records)
        Kp, dna = parse_jellyfish_dna_dump(path)
        @test Kp == 31
        for (seq, count) in records
          bits = Kmer{DNAAlphabet{2}, 31}(LongDNA{2}(seq)).data[1]
          @test dna[bits] == UInt32(count)
        end
      finally
        isfile(path) && rm(path)
      end
    end

    @testset "min_count filtering" begin
      records = [("ATGCCC", 4), ("AAACCC", 200)]
      path = tempname()
      try
        _write_jf_dump(path, K, count_bytes, records)
        _, dna = parse_jellyfish_dna_dump(path; min_count=100)
        @test length(dna) == 1
        @test only(values(dna)) == UInt32(200)
      finally
        isfile(path) && rm(path)
      end
    end

    @testset "canonical dumps rejected unless allow_canonical=true" begin
      path = tempname()
      try
        _write_jf_dump(path, K, count_bytes, []; canonical=true)
        @test_throws ErrorException parse_jellyfish_dna_dump(path)
        _, dna = parse_jellyfish_dna_dump(path; allow_canonical=true)
        @test isempty(dna)
      finally
        isfile(path) && rm(path)
      end
    end

    @testset "K > 32 rejected" begin
      path = tempname()
      try
        _write_jf_dump(path, 33, 4, [])
        @test_throws ErrorException parse_jellyfish_dna_dump(path)
      finally
        isfile(path) && rm(path)
      end
    end

    @testset "jellyfish_dump_hash: translation, stop-codon drop, synonymous-codon summation" begin
      # AAATAA: AAA=Lys, TAA=Stop -> dropped. GCT/GCC both -> Ala, paired with ATG=Met -> summed.
      records = [("ATGCCC", 4), ("AAATAA", 9), ("ATGGCT", 3), ("ATGGCC", 4)]
      path = tempname()
      try
        _write_jf_dump(path, K, count_bytes, records)
        aa = jellyfish_dump_hash(path)
        @test length(aa) == 2
        @test aa[Kmer{AAAlphabet, 2}(LongAA("MP")).data[1]] == UInt32(4)
        @test aa[Kmer{AAAlphabet, 2}(LongAA("MA")).data[1]] == UInt32(7)

        kct = KCT{K÷3, AAAlphabet}(aa)
        @test length(kct) == 2
      finally
        isfile(path) && rm(path)
      end
    end
  end

  # Full count vector for k-mer `x` across `samples`, padded with zeros for absences.
  # This is exactly what kct.counts[pos] must return once all samples are in the table.
  _expected_cv(samples, x) = UInt32[get(s, x, UInt32(0)) for s in samples]

  @testset "8. Batched push! == sequential single push!" begin
    s1 = Dict{UInt64, UInt32}(0x11 => 3, 0x22 => 5, 0x33 => 1)
    s2 = Dict{UInt64, UInt32}(0x11 => 4, 0x44 => 2, 0x55 => 9)
    s3 = Dict{UInt64, UInt32}(0x11 => 2, 0x22 => 7, 0x44 => 6, 0x66 => 8)
    s4 = Dict{UInt64, UInt32}(0x11 => 1, 0x66 => 5, 0x77 => 4)
    samples = (s1, s2, s3, s4)
    allk = sort(collect(UInt64, union(keys.(samples)...)))

    # Sequential: one push! per sample, collapse! between (mirrors the build drivers).
    seq = KCT{4, AAAlphabet}(deepcopy(s1))
    for s in (s2, s3, s4)
      push!(seq, deepcopy(s))
      seq = collapse!(seq)
    end

    # Batched: sample 1 builds the table, the rest go in as one batch.
    bat = KCT{4, AAAlphabet}(deepcopy(s1))
    push!(bat, [deepcopy(s2), deepcopy(s3), deepcopy(s4)])
    bat = collapse!(bat)

    @test bat.counts.samples[] == 4
    @test seq.counts.samples[] == 4
    @test length(bat) == length(allk) == 7
    @test collect(bat.kmer.seqs) == collect(seq.kmer.seqs) == allk

    for x in allk
      want = _expected_cv(samples, x)
      @test bat.counts[findfirst(bat.kmer, x)] == want
      @test seq.counts[findfirst(seq.kmer, x)] == want
    end
  end

  @testset "9. Batched push! edge cases (gaps, mid-batch new, last-sample new)" begin
    # base has one k-mer. The batch exercises: embedded zero (0xA1 in b1 & b3, not b2),
    # trailing implicit zero (0xA2 in b1 only), leading backfill (0xB0 in b2 only),
    # brand-new mid-batch then seen again (0xC0 in b2 & b3), brand-new in the last
    # batch sample only (0xD0 in b3).
    base = Dict{UInt64, UInt32}(0xA1 => 2, 0xA2 => 9)
    b1 = Dict{UInt64, UInt32}(0xA1 => 5)
    b2 = Dict{UInt64, UInt32}(0xB0 => 4, 0xC0 => 7)
    b3 = Dict{UInt64, UInt32}(0xA1 => 6, 0xC0 => 1, 0xD0 => 8)
    all_samples = (base, b1, b2, b3)

    kct = KCT{4, AAAlphabet}(deepcopy(base))
    push!(kct, [deepcopy(b1), deepcopy(b2), deepcopy(b3)])
    kct = collapse!(kct)

    @test kct.counts.samples[] == 4
    allk = sort(collect(UInt64, union(keys.(all_samples)...)))
    @test collect(kct.kmer.seqs) == allk

    for x in allk
      @test kct.counts[findfirst(kct.kmer, x)] == _expected_cv(all_samples, x)
    end

    # spot-check the tricky vectors explicitly
    @test kct.counts[findfirst(kct.kmer, UInt64(0xA1))] == UInt32[2, 5, 0, 6]
    @test kct.counts[findfirst(kct.kmer, UInt64(0xA2))] == UInt32[9, 0, 0, 0]
    @test kct.counts[findfirst(kct.kmer, UInt64(0xB0))] == UInt32[0, 0, 4, 0]
    @test kct.counts[findfirst(kct.kmer, UInt64(0xC0))] == UInt32[0, 0, 7, 1]
    @test kct.counts[findfirst(kct.kmer, UInt64(0xD0))] == UInt32[0, 0, 0, 8]
  end

  @testset "10. Streaming builder (V4.0 SparseCountsLayer)" begin
    # K=12 nucleotides -> 4 AA symbols -> 20-bit encoding. shard_bits=2 (keep_shift=18)
    # spreads these k-mers across all 4 prefix shards.
    K1 = UInt64(0x00100); K2 = UInt64(0x02000); K3 = UInt64(0x40001)
    K4 = UInt64(0x80005); K5 = UInt64(0xC0002); K6 = UInt64(0x40ABC)
    S = [
      Dict{UInt64, UInt32}(K1 => 3, K2 => 5, K3 => 1),
      Dict{UInt64, UInt32}(K1 => 3, K4 => 2, K3 => 1),
      Dict{UInt64, UInt32}(K1 => 2, K2 => 7, K4 => 6, K5 => 8),
      Dict{UInt64, UInt32}(K3 => 4, K6 => 9),
      Dict{UInt64, UInt32}(K5 => 5, K6 => 9),
      Dict{UInt64, UInt32}(K1 => 1, K2 => 5),
    ]
    allk = sort([K1, K2, K3, K4, K5, K6])
    paths = ["s$i" for i in eachindex(S)]
    cmap = Dict(paths[i] => S[i] for i in eachindex(S))
    mkcounter(m) = (p -> deepcopy(m[p]))

    stream_build(; kw...) = mktempdir() do tmp
      mktempdir() do out
        outp = build_kct_streaming(paths; K = 12, translate = true, tmp_dir = tmp, out_dir = out,
                                   counter = mkcounter(cmap), kw...)
        (get_version(outp), load_kct(outp))  # dir is torn down after load
      end
    end

    @testset "streaming == incremental; hierarchical merge" begin
      # wave_size=2 -> 3 waves, merge_fanin=2 -> a real merge level before assembly.
      ver, kct = stream_build(wave_size = 2, shard_bits = 2, merge_fanin = 2)
      @test ver == 4.0
      @test kct.counts isa SparseCountsLayer
      @test length(kct) == 6
      @test kct.counts.n_samples.x == 6
      @test collect(kct.kmer.seqs) == allk
      for x in allk
        @test kct.counts[findfirst(kct.kmer, x)] == _expected_cv(S, x)
      end
      # K3 spans non-adjacent waves 1 and 2 (samples 1,2,4), K6 spans waves 2,3 (4,5)
      @test kct.counts[findfirst(kct.kmer, K3)] == UInt32[1, 1, 0, 4, 0, 0]
      @test kct.counts[findfirst(kct.kmer, K6)] == UInt32[0, 0, 0, 9, 9, 0]
    end

    @testset "single-wave path (no merge levels)" begin
      ver, kct = stream_build(wave_size = 100, shard_bits = 2, merge_fanin = 5)
      @test ver == 4.0
      for x in allk
        @test kct.counts[findfirst(kct.kmer, x)] == _expected_cv(S, x)
      end
    end

    @testset "sparsity: only present samples are stored" begin
      _, kct = stream_build(wave_size = 2, shard_bits = 2, merge_fanin = 2)
      scl = kct.counts
      # K5 is present in exactly samples 3 and 5 -> its row has exactly 2 pairs
      r = scl.row_id[findfirst(kct.kmer, K5)]
      @test scl.row_offsets[r + 1] - scl.row_offsets[r] == 2
      # K1 present in 4 of 6 samples -> 4 pairs (no zeros stored)
      r1 = scl.row_id[findfirst(kct.kmer, K1)]
      @test scl.row_offsets[r1 + 1] - scl.row_offsets[r1] == 4
    end

    @testset "row deduplication" begin
      # A,B,C,D all == count 4 in sample 1 only, so they must collapse to one pooled row.
      A = UInt64(0x00001); B = UInt64(0x00002); C = UInt64(0x40001)
      Dk = UInt64(0x80001); E = UInt64(0xC0001)
      d1 = Dict{UInt64, UInt32}(A => 4, B => 4, C => 4, Dk => 4, E => 9)
      d2 = Dict{UInt64, UInt32}(E => 9)
      cm = Dict("a" => d1, "b" => d2)
      kct = mktempdir() do tmp; mktempdir() do out
        load_kct(build_kct_streaming(["a", "b"]; K = 12, translate = true, wave_size = 1, shard_bits = 2,
                                     merge_fanin = 2, tmp_dir = tmp, out_dir = out,
                                     counter = (p -> deepcopy(cm[p]))))
      end; end
      rid = x -> kct.counts.row_id[findfirst(kct.kmer, x)]
      @test rid(A) == rid(B) == rid(C) == rid(Dk)
      @test rid(E) != rid(A)
      @test length(kct.counts.row_offsets) - 1 == 2  # exactly two distinct rows
      for x in (A, B, C, Dk, E)
        @test kct.counts[findfirst(kct.kmer, x)] == UInt32[get(d1, x, UInt32(0)), get(d2, x, UInt32(0))]
      end
    end

    @testset "V3.0 seed adapter" begin
      full = [
        Dict{UInt64, UInt32}(K1 => 3, K2 => 5),
        Dict{UInt64, UInt32}(K1 => 2, K3 => 7),
        Dict{UInt64, UInt32}(K2 => 1, K4 => 9),
        Dict{UInt64, UInt32}(K1 => 4, K3 => 2, K4 => 1),
      ]
      seed = KCT{4, AAAlphabet}(deepcopy(full[1]))
      push!(seed, deepcopy(full[2])); seed = collapse!(seed)

      kct = mktempdir() do sd; mktempdir() do tmp; mktempdir() do out
        sp = joinpath(sd, "seed.kct"); write_kct(seed, sp)
        cm = Dict("s3" => full[3], "s4" => full[4])
        load_kct(build_kct_streaming(["s3", "s4"]; K = 12, translate = true, wave_size = 1, shard_bits = 2,
                                     merge_fanin = 2, tmp_dir = tmp, out_dir = out,
                                     counter = (p -> deepcopy(cm[p])), seeds = [(sp, 1)]))
      end; end; end
      @test kct.counts.n_samples.x == 4
      @test length(kct) == 4
      for x in sort([K1, K2, K3, K4])
        @test kct.counts[findfirst(kct.kmer, x)] == UInt32[get(s, x, UInt32(0)) for s in full]
      end
    end
  end

  @testset "11. Rolling k-mer counter and DNA streaming build" begin
    _count(chunk, K; translate) = begin
      mq = Channel{Dict{UInt64, UInt32}}(1)
      count_kmers(chunk, K, mq; translate = translate)
      take!(mq)
    end
    _ok(l) = length(l) >= 12 && all(c -> c in ('A', 'C', 'G', 'T'), l)

    reads = ["ACGTACGTACGTACGTACGTACGTACGTAC",
             "GGGCCCTTTAAAGGGCCCTTTAAAGGGCCCT",
             "ACGTNCGTACGTACGTACGTACGTACGTAC",   # N -> whole read skipped
             "ACG",                              # shorter than K -> skipped
             "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
             "ATGATGATGATGATGATGATGATGATGATG"]
    K = 12

    @testset "translate=true is bit-identical to translate(k_merize(...))" begin
      ref = Dict{UInt64, UInt32}()
      for l in filter(_ok, reads)
        for km in k_merize(LongSequence{DNAAlphabet{2}}(l), K = K)
          a = translate(km)
          a === nothing && continue
          ref[a.data[1]] = get(ref, a.data[1], UInt32(0)) + UInt32(1)
        end
      end
      @test _count(reads, K; translate = true) == ref
    end

    @testset "translate=false gives raw DNA k-mer codes" begin
      ref = Dict{UInt64, UInt32}()
      for l in filter(_ok, reads)
        seq = LongSequence{DNAAlphabet{2}}(l)
        for i in 1:(length(seq) - K + 1)
          c = Kmer{DNAAlphabet{2}, K}(seq[i:i + K - 1]).data[1]
          ref[c] = get(ref, c, UInt32(0)) + UInt32(1)
        end
      end
      @test _count(reads, K; translate = false) == ref
    end

    @testset "DNA streaming build -> KCT{K, DNAAlphabet{2}} V4.0" begin
      A = UInt64(0x000010); B = UInt64(0x200001); C = UInt64(0x400005)
      Dk = UInt64(0x600002); E = UInt64(0x400777)
      S = [
        Dict{UInt64, UInt32}(A => 3, B => 5, C => 1),
        Dict{UInt64, UInt32}(A => 3, Dk => 2, C => 1),
        Dict{UInt64, UInt32}(A => 2, B => 7, Dk => 6, E => 8),
        Dict{UInt64, UInt32}(C => 4, E => 9),
      ]
      allk = sort([A, B, C, Dk, E])
      cv(x) = UInt32[get(s, x, UInt32(0)) for s in S]
      cm = Dict("s$i" => S[i] for i in eachindex(S))

      ver, kct = mktempdir() do tmp; mktempdir() do out
        p = build_kct_streaming(["s$i" for i in eachindex(S)]; K = 12, translate = false, idx_prefix = 12,
                                wave_size = 2, shard_bits = 2, merge_fanin = 2,
                                tmp_dir = tmp, out_dir = out, counter = (x -> deepcopy(cm[x])))
        (get_version(p), load_kct(p))
      end; end

      @test ver == 4.0
      @test kct isa KCT{12, DNAAlphabet{2}}
      @test kct.counts isa SparseCountsLayer
      @test kct.counts.n_samples.x == 4
      @test collect(kct.kmer.seqs) == allk
      for x in allk
        @test kct.counts[findfirst(kct.kmer, x)] == cv(x)
      end

      # round-trip preserves the DNA alphabet in the header
      kct2 = mktempdir() do d
        p = joinpath(d, "rt.kct"); write_kct(kct, p); load_kct(p)
      end
      @test kct2 isa KCT{12, DNAAlphabet{2}}
      for x in allk
        @test kct2.counts[findfirst(kct2.kmer, x)] == cv(x)
      end
    end
  end

end
