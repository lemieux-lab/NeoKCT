using Pkg
Pkg.activate(".")

include("NeoKCT.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    samples = readlines(open("/u/jacquinn/phd_stuff/data/tcga_fastqs.paths", "r"))
    build_kct(samples, 
    save_at_samples= [5, 10, 20, 50, 100, 200, 500, 750, 1000, collect(1500:500:10935)..., 10935], 
    save_path="/u/jacquinn/phd_stuff/data/neo_kcts/production/TCGA_V1.0/",
    word_size=UInt64)
end