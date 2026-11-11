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

# Full count vector for k-mer `x` across `samples`, padded with zeros for absences.
# This is exactly what kct.counts[pos] must return once all samples are in the table.
_expected_cv(samples, x) = UInt32[get(s, x, UInt32(0)) for s in samples]

@testset "NeoKCT" begin

  @testset "1. BiotypLayer bitmask logic" begin
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

  @testset "2. KCTLoader write/load round-trip" begin
    raw = Dict{UInt64, UInt32}(
      UInt64(100) => UInt32(3),
      UInt64(200) => UInt32(7),
      UInt64(500) => UInt32(1),
    )
    kct = mktempdir() do tmp; mktempdir() do out
      p = build_kct_streaming(["s1"]; K = 12, translate = true, tmp_dir = tmp, out_dir = out,
                              counter = (_ -> deepcopy(raw)))
      load_kct(p)
    end; end

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

  @testset "3. JellyfishDump parsing" begin
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
      finally
        isfile(path) && rm(path)
      end
    end
  end

  @testset "4. Streaming builder (V4.0 CountsLayer)" begin
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
      @test kct.counts isa CountsLayer
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
      # K5 is present in exactly samples 3 and 5 -> its reconstructed vector has 2 non-zeros
      v5 = kct.counts[findfirst(kct.kmer, K5)]
      @test length(v5) == 6
      @test count(!iszero, v5) == 2
      # K1 present in 4 of 6 samples -> 4 non-zeros, no stored zeros
      @test count(!iszero, kct.counts[findfirst(kct.kmer, K1)]) == 4
    end

    @testset "identical count vectors reconstruct identically" begin
      # No cross-k-mer dedup; identical vectors must still round-trip identically.
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
      cv = x -> kct.counts[findfirst(kct.kmer, x)]
      @test cv(A) == cv(B) == cv(C) == cv(Dk) == UInt32[4, 0]
      @test cv(E) == UInt32[9, 9]
      for x in (A, B, C, Dk, E)
        @test cv(x) == UInt32[get(d1, x, UInt32(0)), get(d2, x, UInt32(0))]
      end
    end

    @testset "seed adapter" begin
      full = [
        Dict{UInt64, UInt32}(K1 => 3, K2 => 5),
        Dict{UInt64, UInt32}(K1 => 2, K3 => 7),
        Dict{UInt64, UInt32}(K2 => 1, K4 => 9),
        Dict{UInt64, UInt32}(K1 => 4, K3 => 2, K4 => 1),
      ]
      kct = mktempdir() do seed_tmp; mktempdir() do seed_out; mktempdir() do tmp; mktempdir() do out
        seed_cm = Dict("s1" => full[1], "s2" => full[2])
        sp = build_kct_streaming(["s1", "s2"]; K = 12, translate = true, wave_size = 2, shard_bits = 2,
                                 merge_fanin = 2, tmp_dir = seed_tmp, out_dir = seed_out,
                                 counter = (p -> deepcopy(seed_cm[p])))
        cm = Dict("s3" => full[3], "s4" => full[4])
        load_kct(build_kct_streaming(["s3", "s4"]; K = 12, translate = true, wave_size = 1, shard_bits = 2,
                                     merge_fanin = 2, tmp_dir = tmp, out_dir = out,
                                     counter = (p -> deepcopy(cm[p])), seeds = [(sp, 1)]))
      end; end; end; end
      @test kct.counts.n_samples.x == 4
      @test length(kct) == 4
      for x in sort([K1, K2, K3, K4])
        @test kct.counts[findfirst(kct.kmer, x)] == UInt32[get(s, x, UInt32(0)) for s in full]
      end
    end

    @testset "adaptive reshard" begin
      # 90 bytes is verified (not guessed) against this fixture's real shard sizes: every
      # base shard's combined-across-2-roots size is 100-200 bytes, safely above it, so every
      # shard actually splits at least once; an empty child's combined floor is 64 bytes,
      # safely below it, so empty branches terminate in one check instead of recursing to
      # the bit-exhaustion depth limit.
      @testset "forced split (default fanout) reproduces the un-resharded result" begin
        ver, kct = stream_build(wave_size = 2, shard_bits = 2, merge_fanin = 2,
                                reshard_leaf_bytes = 90)
        @test ver == 4.0
        @test collect(kct.kmer.seqs) == allk
        for x in allk
          @test kct.counts[findfirst(kct.kmer, x)] == _expected_cv(S, x)
        end
      end

      @testset "deep recursion with a small fanout still reproduces the result" begin
        # K3=0x40001 and K6=0x40ABC share an 8-bit prefix beyond shard_bits=2, so a 2-way
        # fanout needs several levels to separate them -- exercises real recursion depth,
        # not just one split.
        ver, kct = stream_build(wave_size = 2, shard_bits = 2, merge_fanin = 2,
                                reshard_leaf_bytes = 90, reshard_extra_bits = 1)
        @test ver == 4.0
        @test collect(kct.kmer.seqs) == allk
        for x in allk
          @test kct.counts[findfirst(kct.kmer, x)] == _expected_cv(S, x)
        end
      end

      @testset "_split_shard_file! preserves every (kmer, sample, count) triple" begin
        mktempdir() do dir
          # kLow and kHigh share a pooled row (identical sample/count pattern) despite
          # landing in DIFFERENT children -- the real-world case a naive memoized
          # old_row_id -> new_row_id shortcut would misassign.
          rp = _RowPool()
          shared = Tuple{UInt32, UInt32}[(UInt32(1), UInt32(9))]
          r_shared = _intern_row!(rp, shared)
          pairs_mid = Tuple{UInt32, UInt32}[(UInt32(2), UInt32(3)), (UInt32(5), UInt32(7))]
          r_mid = _intern_row!(rp, pairs_mid)

          kLow = UInt64(1)    # top bit (of a 5-bit window) 0 -> child 0
          kMid = UInt64(15)   # top bit 0 -> child 0, own row
          kHigh = UInt64(17)  # top bit 1 -> child 1, shares kLow's row
          ps = _PartialShard(UInt64[kLow, kMid, kHigh], UInt32[r_shared, r_mid, r_shared],
                             rp.offsets, rp.pool_s, rp.pool_c)
          src = joinpath(dir, "src.bin")
          _write_partial_shard(src, ps)

          _split_shard_file!(src, dir, "x", 5, 1)

          expected = Set{Tuple{UInt64, UInt32, UInt32}}()
          for p in shared; push!(expected, (kLow, p[1], p[2])); end
          for p in pairs_mid; push!(expected, (kMid, p[1], p[2])); end
          for p in shared; push!(expected, (kHigh, p[1], p[2])); end

          got = Set{Tuple{UInt64, UInt32, UInt32}}()
          for c in 0:1
            cps = _read_partial_shard(joinpath(dir, "x_$(c).bin"))
            for i in eachindex(cps.kmers)
              r = cps.row_id[i]
              for t in cps.offsets[r]:(cps.offsets[r + 1] - 1)
                push!(got, (cps.kmers[i], cps.pool_s[t], cps.pool_c[t]))
              end
            end
          end
          @test got == expected
        end
      end

      @testset "resumability: redoing a prefix doesn't mix stale children" begin
        mktempdir() do outer; mktempdir() do rootdir
          ps = _PartialShard(UInt64[0b00, 0b01, 0b10, 0b11], UInt32[1, 2, 3, 4],
                             UInt64[1, 2, 3, 4, 5], UInt32[10, 20, 30, 40], UInt32[1, 1, 1, 1])
          _write_partial_shard(joinpath(rootdir, "shard0.bin"), ps)
          work = joinpath(outer, "work"); mkpath(work)
          mark_dir = joinpath(work, "_markers"); mkpath(mark_dir)
          _done(name) = isfile(joinpath(mark_dir, name))
          _mark(name) = touch(joinpath(mark_dir, name))

          leaves1 = _reshard_leaves([rootdir], 1, 2, work, 1, 1, _done, _mark)
          @test length(leaves1) == 4  # keep_shift=2, extra_bits=1 -> exactly 2 levels

          # Simulate a crash recovered after the top split wrote its children but before it
          # was marked: redo must only re-split this prefix, not disturb already-marked
          # deeper levels, and must still land on the same final leaf set.
          rm(joinpath(mark_dir, "reshard_shard0"))
          leaves2 = _reshard_leaves([rootdir], 1, 2, work, 1, 1, _done, _mark)
          @test leaves2 == leaves1
          for leafpaths in leaves2
            @test all(isfile, leafpaths)
          end
        end; end
      end

      @testset "reshard_extra_bits is pinned across resumed runs" begin
        mktempdir() do tmp; mktempdir() do out
          # A huge leaf target means nothing actually splits here -- this call only needs to
          # reach the reshard stage once successfully, so the pin file gets written with
          # reshard_extra_bits=4 (a tiny leaf_bytes, as used elsewhere in this testset to
          # force real splitting, would be wasteful and pointless for this test: the pin
          # file is written unconditionally before any splitting decision is made).
          build_kct_streaming(["s1"]; K = 12, translate = true, wave_size = 1, shard_bits = 2,
                              tmp_dir = tmp, out_dir = out, counter = (p -> deepcopy(S[1])),
                              reshard_leaf_bytes = 1_000_000_000, reshard_extra_bits = 4,
                              resume = true)
          # A fresh call against the SAME out_dir with a different reshard_extra_bits must
          # error rather than silently misreading the bit position already-split children
          # were written at (work was deleted by the successful run above, so re-seed it
          # with one marker to simulate an interrupted run that got as far as pinning it).
          work = joinpath(out, "streaming_work"); mark_dir = joinpath(work, "_markers")
          mkpath(mark_dir)
          write(joinpath(mark_dir, "reshard_extra_bits.n"), "4")
          @test_throws ErrorException build_kct_streaming(
              ["s1"]; K = 12, translate = true, wave_size = 1, shard_bits = 2,
              tmp_dir = tmp, out_dir = out, counter = (p -> deepcopy(S[1])),
              reshard_leaf_bytes = 1, reshard_extra_bits = 1, resume = true)
        end; end
      end
    end

    @testset "patch_sample_counts" begin
      K_new = UInt64(0x40500)  # between K3 and K6, present in neither S nor allk
      corrected = Dict{UInt64, UInt32}(K1 => 99, K3 => 1, K_new => 42)  # no K2 -> scrub test

      kct2 = mktempdir() do tmp; mktempdir() do out; mktempdir() do patchdir
        src = build_kct_streaming(paths; K = 12, translate = true, tmp_dir = tmp, out_dir = out,
                                  counter = mkcounter(cmap), wave_size = 2, shard_bits = 2, merge_fanin = 2)
        dst = joinpath(patchdir, "patched.kct")
        patch_sample_counts(src, 1, corrected, dst)
        load_kct(dst)
      end; end; end

      @test kct2.counts.n_samples.x == 6  # patching never adds a sample slot
      @test collect(kct2.kmer.seqs) == sort(vcat(allk, [K_new]))

      expected1 = copy(_expected_cv(S, K1)); expected1[1] = 99  # touched: count replaced
      @test kct2.counts[findfirst(kct2.kmer, K1)] == expected1

      expected2 = copy(_expected_cv(S, K2)); expected2[1] = 0  # scrubbed: absent from corrected
      @test kct2.counts[findfirst(kct2.kmer, K2)] == expected2

      @test kct2.counts[findfirst(kct2.kmer, K3)] == _expected_cv(S, K3)  # same value either way

      expectednew = zeros(UInt32, 6); expectednew[1] = 42  # brand-new k-mer, one-sample row
      @test kct2.counts[findfirst(kct2.kmer, K_new)] == expectednew

      for x in (K4, K5, K6)  # sample 1 never touched these -- must round-trip unchanged
        @test kct2.counts[findfirst(kct2.kmer, x)] == _expected_cv(S, x)
      end

      @testset "a row scrubbed down to zero samples is dropped entirely" begin
        Ksolo = UInt64(0x00300)  # present only in sample "a", absent from the corrected recount
        Kboth = UInt64(0x00400)  # present in both, untouched by the patch
        cm2 = Dict("a" => Dict{UInt64, UInt32}(Ksolo => 7, Kboth => 2),
                   "b" => Dict{UInt64, UInt32}(Kboth => 3))
        kct3 = mktempdir() do tmp; mktempdir() do out; mktempdir() do patchdir
          src = build_kct_streaming(["a", "b"]; K = 12, translate = true, wave_size = 1, shard_bits = 2,
                                    merge_fanin = 2, tmp_dir = tmp, out_dir = out,
                                    counter = (p -> deepcopy(cm2[p])))
          dst = joinpath(patchdir, "patched.kct")
          patch_sample_counts(src, 1, Dict{UInt64, UInt32}(Kboth => 2), dst)
          load_kct(dst)
        end; end; end
        @test collect(kct3.kmer.seqs) == [Kboth]
        @test kct3.counts[findfirst(kct3.kmer, Kboth)] == UInt32[2, 3]
      end

      @testset "out-of-range sample_slot errors" begin
        mktempdir() do tmp; mktempdir() do out; mktempdir() do patchdir
          src = build_kct_streaming(paths; K = 12, translate = true, tmp_dir = tmp, out_dir = out,
                                    counter = mkcounter(cmap), wave_size = 2, shard_bits = 2, merge_fanin = 2)
          dst = joinpath(patchdir, "patched.kct")
          @test_throws ErrorException patch_sample_counts(src, 0, Dict{UInt64, UInt32}(), dst)
          @test_throws ErrorException patch_sample_counts(src, 7, Dict{UInt64, UInt32}(), dst)
        end; end; end
      end
    end
  end

  @testset "5. Rolling k-mer counter and DNA streaming build" begin
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

    @testset "double_strand=false (default) is unchanged" begin
      mq = Channel{Dict{UInt64, UInt32}}(1)
      count_kmers(reads, K, mq; translate = false)
      @test take!(mq) == _count(reads, K; translate = false)
    end

    @testset "double_strand=true counts each read's reverse complement too" begin
      _count_ds(chunk, K; translate) = begin
        mq = Channel{Dict{UInt64, UInt32}}(1)
        count_kmers(chunk, K, mq; translate = translate, double_strand = true)
        take!(mq)
      end
      _revcomp_ref(l) = string(BioSequences.reverse_complement(LongDNA{4}(l)))

      @testset "translate=false" begin
        forward = _count(reads, K; translate = false)
        revcomp = _count(map(_revcomp_ref, filter(_ok, reads)), K; translate = false)
        @test _count_ds(reads, K; translate = false) == merge(+, forward, revcomp)
      end

      @testset "translate=true" begin
        forward = _count(reads, K; translate = true)
        revcomp = _count(map(_revcomp_ref, filter(_ok, reads)), K; translate = true)
        @test _count_ds(reads, K; translate = true) == merge(+, forward, revcomp)
      end

      @testset "a palindromic-adjacent read still counts both strands independently" begin
        # ACGT's own reverse complement is ACGT, so a read built entirely from it is its own
        # reverse complement -- double_strand should still double its counts, not collapse
        # them, since count_strand! is called twice regardless of what the two strands equal.
        pal = ["ACGTACGTACGT"]
        K2 = 4
        once = _count(pal, K2; translate = false)
        @test _count_ds(pal, K2; translate = false) == Dict(k => 2v for (k, v) in once)
      end
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
      @test kct.counts isa CountsLayer
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

  @testset "12. Peptide queries" begin
    # K=18nt -> 6aa keeps DEFAULT_IDX_PREFIX_SIZE (5) meaningfully smaller than K, so the
    # prefix index partitions on exactly the first residue (32 buckets, bucket = top 5 bits
    # = 1 symbol) instead of collapsing to a single bucket the way the K=4aa tables above do.
    # Real bucket boundaries, not just correct arithmetic, is the point of this testset.
    kk = 6

    # M-bucket: 4 k-mers sharing first residue M, sorted "MA.." before "MC..".
    K1 = "MAAAAA"; K2 = "MAACCC"; K3 = "MCCCCC"; K4 = "MCCCCG"
    # K-bucket: 2 k-mers, a different first residue, isolation check against the M-bucket.
    K5 = "KAAAAA"; K6 = "KCCCCC"
    # Second tile for an L > K peptide ("MAAAAA" + "AAAAAG" overlap by 5 residues).
    K7 = "AAAAAG"
    # I/L pair at the same position, same counts pattern as K1/K2 in shape but distinct rows.
    K8 = "MLAAAA"; K9 = "MIAAAA"

    peps = [K1, K2, K3, K4, K5, K6, K7, K8, K9]
    counts = Dict(
      K1 => (3, 0, 5), K2 => (2, 4, 0), K3 => (1, 1, 1), K4 => (0, 2, 3),
      K5 => (7, 0, 0), K6 => (0, 6, 0), K7 => (4, 4, 4), K8 => (2, 3, 1), K9 => (5, 0, 2),
    )
    codes = Dict(p => _aa_code(p) for p in peps)
    S = [Dict{UInt64, UInt32}(codes[p] => counts[p][s] for p in peps if counts[p][s] != 0) for s in 1:3]
    paths = ["s1", "s2", "s3"]
    cmap = Dict(paths[i] => S[i] for i in eachindex(S))

    kct = mktempdir() do tmp; mktempdir() do out
      p = build_kct_streaming(paths; K = 3kk, translate = true, wave_size = 3, shard_bits = 2,
                              tmp_dir = tmp, out_dir = out, counter = (x -> deepcopy(cmap[x])))
      load_kct(p)
    end; end
    @test idx_prefix_size(kct.kmer) == 5

    @testset "_aa_code matches translate() encoding" begin
      dna = LongSequence{DNAAlphabet{2}}("ATGGCGGCGGCGGCAGCC")  # 18nt, 6 codons M-A-A-A-A-A, no stop
      km = Kmer{DNAAlphabet{2}, 18}(dna)
      aa_km = translate(km)
      @test aa_km !== nothing
      pep_str = String([Char(aa_km[i]) for i in 1:kk])
      @test _aa_code(pep_str) == aa_km.data[1]
    end

    @testset "L < K: prefix range sums, real bucket boundaries" begin
      # "M" spans the whole M-bucket: first_hit lands at the bucket's own r.start.
      # (K1+K2+K3+K4+K8+K9 -- the I/L pair also starts with M, same bucket as K1-K4.)
      @test peptide_counts(kct, "M") == UInt32[13, 10, 12]
      # "MA" is the bucket's leading sub-run (K1, K2).
      @test peptide_counts(kct, "MA") == UInt32[5, 4, 5]
      # "MC" is the bucket's trailing sub-run (K3, K4), ending exactly at r.stop.
      @test peptide_counts(kct, "MC") == UInt32[1, 3, 4]
      # "K" is a different bucket entirely: isolation from the M-bucket's contents.
      @test peptide_counts(kct, "K") == UInt32[7, 6, 0]
    end

    @testset "L == K: exact match" begin
      @test peptide_counts(kct, K1) == UInt32[3, 0, 5]
    end

    @testset "L > K: tile minimum" begin
      # "MAAAAAG" tiles as K1="MAAAAA" + K7="AAAAAG"; bounds the count from above.
      @test peptide_counts(kct, "MAAAAAG") == min.(UInt32[3, 0, 5], UInt32[4, 4, 4])
    end

    @testset "C-terminal blind spot: absent prefix returns zero, not an error" begin
      # No k-mer in this table starts with residue "G"; documents the blind spot rather
      # than fixing it (PAPER_TODO.md task A5).
      @test peptide_counts(kct, "G") == UInt32[0, 0, 0]
    end

    @testset "il_ambiguous sums every I/L variant" begin
      @test peptide_counts(kct, "MLAAAA") == UInt32[2, 3, 1]
      @test peptide_counts(kct, "MIAAAA") == UInt32[5, 0, 2]
      @test peptide_counts(kct, "MLAAAA"; il_ambiguous = true) == UInt32[7, 3, 3]
      @test peptide_counts(kct, "MIAAAA"; il_ambiguous = true) == UInt32[7, 3, 3]
    end

    @testset "peptide_matrix matches per-peptide peptide_counts" begin
      qpeps = [K1, "M", "MA", "MAAAAAG", "G"]
      M = peptide_matrix(kct, qpeps)
      for (j, p) in enumerate(qpeps)
        @test M[j, :] == peptide_counts(kct, p)
      end
    end

    @testset "searchsorted range-spanning-buckets throws rather than misreading" begin
      # idx_prefix_size(kct.kmer) == 5 for this table, so a 0-residue "prefix" (free = 30
      # bits) cannot fit in one bucket: exercises the ArgumentError guard directly.
      @test_throws ArgumentError searchsorted(kct.kmer, UInt64(0), typemax(UInt64) >> (64 - 30))
    end
  end

  @testset "13. Joins (_merge_walk, add_biotypes, setdiff)" begin
    @testset "_merge_walk: exact (i, j) pairs on known overlapping arrays" begin
      a = DeltaArray(UInt64[1, 3, 5, 7, 9])
      b = DeltaArray(UInt64[3, 5, 6, 9, 10])
      pairs = Tuple{Int, Int}[]
      _merge_walk((i, j) -> push!(pairs, (i, j)), a, b)
      # value 3: a[2], b[1]. value 5: a[3], b[2]. value 9: a[5], b[4].
      @test pairs == [(2, 1), (3, 2), (5, 4)]
    end

    @testset "add_biotypes: matched k-mers get the right mask, unmatched stay intergenic" begin
      bnames = ["intergenic", "protein_coding", "lncRNA"]
      gidx_kmers = UInt64[0x10, 0x20, 0x30, 0x40]
      gidx_masks = UInt64[2, 4, 2, 4]  # 0x10,0x30 -> protein_coding; 0x20,0x40 -> lncRNA
      gidx = KCT{4, AAAlphabet}(gidx_kmers, gidx_masks, bnames)

      S = [Dict{UInt64, UInt32}(UInt64(0x10) => 3, UInt64(0x25) => 5),
           Dict{UInt64, UInt32}(UInt64(0x30) => 7, UInt64(0x50) => 2)]
      paths = ["s1", "s2"]
      cmap = Dict(paths[i] => S[i] for i in eachindex(S))
      kct = mktempdir() do tmp; mktempdir() do out
        p = build_kct_streaming(paths; K = 12, translate = true, wave_size = 2, shard_bits = 1,
                                tmp_dir = tmp, out_dir = out, counter = (x -> deepcopy(cmap[x])))
        load_kct(p)
      end; end
      @test collect(kct.kmer.seqs) == sort(UInt64[0x10, 0x25, 0x30, 0x50])

      rich = add_biotypes(kct, gidx)
      @test has_biotype(rich.biotype, findfirst(rich.kmer, UInt64(0x10)), "protein_coding")
      @test has_biotype(rich.biotype, findfirst(rich.kmer, UInt64(0x30)), "protein_coding")
      @test has_biotype(rich.biotype, findfirst(rich.kmer, UInt64(0x25)), "intergenic")
      @test has_biotype(rich.biotype, findfirst(rich.kmer, UInt64(0x50)), "intergenic")
      # counts pass through untouched
      @test rich.counts[findfirst(rich.kmer, UInt64(0x10))] == kct.counts[findfirst(kct.kmer, UInt64(0x10))]
    end

    @testset "setdiff: default predicate is plain set difference" begin
      tS = [Dict{UInt64, UInt32}(UInt64(0x10) => 4, UInt64(0x20) => 1, UInt64(0x30) => 9)]
      nS = [Dict{UInt64, UInt32}(UInt64(0x20) => 2, UInt64(0x40) => 6)]
      tumor = mktempdir() do tmp; mktempdir() do out
        p = build_kct_streaming(["t1"]; K = 12, translate = true, wave_size = 1, shard_bits = 1,
                                tmp_dir = tmp, out_dir = out, counter = (_ -> deepcopy(tS[1])))
        load_kct(p)
      end; end
      normal = mktempdir() do tmp; mktempdir() do out
        p = build_kct_streaming(["n1"]; K = 12, translate = true, wave_size = 1, shard_bits = 1,
                                tmp_dir = tmp, out_dir = out, counter = (_ -> deepcopy(nS[1])))
        load_kct(p)
      end; end

      idx = setdiff(tumor, normal)
      found = sort([tumor.kmer[i].data[1] for i in idx])
      @test found == sort(UInt64[0x10, 0x30])  # 0x20 is in both, excluded

      # count-aware predicate: keep tumor k-mers with tumor count >= 5 and (absent from
      # normal, or all-zero there) -- 0x30 (tumor=9, absent) passes, 0x10 (tumor=4) does not.
      idx2 = setdiff(tumor, normal; pred = (trow, nrow) -> trow[1] >= 5 && (isnothing(nrow) || all(iszero, nrow)))
      found2 = sort([tumor.kmer[i].data[1] for i in idx2])
      @test found2 == UInt64[0x30]
    end
  end

end
