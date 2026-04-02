using Pkg
Pkg.activate(".")

include("NeoKCT.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    samples = readlines(open("/u/jacquinn/phd_stuff/data/tcga_fastqs.paths", "r"))
    build_kct(samples, 
    benchmark_every = 1,
    save_at_samples= [3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 400, 500, 750, 1000, collect(1500:500:10935)..., 10934, 10935], 
    save_path="/u/jacquinn/phd_stuff/data/neo_kcts/experimental/new_benchmarking_test/",  # Always finish with /
    word_size=UInt64)
end