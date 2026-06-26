using Test
using Kmers, BioSequences, BioSymbols, NArrays, ProgressMeter

include("../src/KCTLayers.jl")

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
      UInt64(100)  => UInt32(3),
      UInt64(200)  => UInt32(7),
      UInt64(500)  => UInt32(1),
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

end
