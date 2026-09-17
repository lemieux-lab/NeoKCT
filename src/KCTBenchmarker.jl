using JSON

# Loading previous history file
function load_benchmark_history(benchmark_file::String)::Vector{Any}
    return isfile(benchmark_file) ? JSON.parsefile(benchmark_file) : []
end

# Adding to / writing history file
function save_benchmark_history!(history::Vector{Any}, entry::Dict{String, Any}, benchmark_file::String)
    push!(history, entry)
    mkpath(dirname(benchmark_file))
    open(benchmark_file, "w") do io
        JSON.print(io, history, 2)
    end
end

"""
    benchmark_kct(kct, benchmark_path; full_pointer_walkthrough=false, benchmark_size=100_000_000)

Measure a KCT's component sizes and k-mer query speed, and append them to the
JSON history under `benchmark_path`. Rows are packed inline block-FOR with no
cross-k-mer dedup, so this reports the k-mer sequence store, the block index,
the packed blob, and query speed. `full_pointer_walkthrough` is accepted for
call-site compatibility with the old V3.0 benchmarker but unused here.
"""
function benchmark_kct(kct::KCT{K, Ab, CountsLayer}, benchmark_path::String;
                       full_pointer_walkthrough::Bool=false, benchmark_size::Int=100_000_000) where {K, Ab<:Alphabet}
    benchmark_file = benchmark_path * "benchmark_data_v4.json"
    history = load_benchmark_history(benchmark_file)
    cl = kct.counts

    printstyled("Measuring KCT components...\n", color=:green)
    n_samples = Int(cl.n_samples.x)
    n_kmers = length(kct.kmer.seqs)
    n_blocks = length(cl.block_ptr)
    checkpoints_bytes = sizeof(eltype(kct.kmer.seqs.checkpoints)) * length(kct.kmer.seqs.checkpoints)
    deltas_bytes = sizeof(eltype(kct.kmer.seqs.deltas)) * length(kct.kmer.seqs.deltas)
    regular_cp_idx_bytes = sizeof(eltype(kct.kmer.seqs.regular_cp_idx)) * length(kct.kmer.seqs.regular_cp_idx)
    kmer_seq_bytes = checkpoints_bytes + deltas_bytes + regular_cp_idx_bytes
    block_ptr_bytes = sizeof(UInt64) * n_blocks
    blob_bytes = length(cl.blob)
    total_bytes = kmer_seq_bytes + block_ptr_bytes + blob_bytes
    bytes_per_kmer = total_bytes / max(n_kmers, 1)

    # count pairs cheaply by decoding a sample of blocks and extrapolating
    sampled = min(n_blocks, 2000)
    pair_sample = 0
    for b in (sampled == 0 ? (1:0) : round.(Int, range(0, n_blocks - 1; length = sampled)))
        for (s, _) in _decode_block(cl, b)
            pair_sample += length(s)
        end
    end
    est_pairs = sampled == 0 ? 0 : round(Int, pair_sample * (n_blocks / sampled))

    k_mers = rand(kct.kmer.seqs.checkpoints, min(benchmark_size, max(1, n_kmers)))
    t_start = now()
    @showprogress "Benchmarking query speed for $(length(k_mers)) k-mers..." for k in k_mers
        findfirst(kct, k)
    end
    query_time_ms = Dates.value(now() - t_start)

    printstyled(
        "$n_kmers k-mers, $n_samples samples, ~$est_pairs pairs, $n_blocks blocks, " *
        "$(round(bytes_per_kmer; digits=2)) B/k-mer\n  kmer_seqs=$(Base.format_bytes(kmer_seq_bytes))  " *
        "block_ptr=$(Base.format_bytes(block_ptr_bytes))  blob=$(Base.format_bytes(blob_bytes))  " *
        "total=$(Base.format_bytes(total_bytes))\n", color=:green)

    entry = Dict{String, Any}(
        "samples" => n_samples, "timestamp" => string(now()),
        "n_kmers" => n_kmers, "n_blocks" => n_blocks, "est_pairs" => est_pairs,
        "kmer_seq_bytes" => kmer_seq_bytes, "block_ptr_bytes" => block_ptr_bytes,
        "blob_bytes" => blob_bytes, "total_bytes" => total_bytes,
        "bytes_per_kmer" => bytes_per_kmer, "query_time_ms" => query_time_ms,
    )
    save_benchmark_history!(history, entry, benchmark_file)
    return
end
