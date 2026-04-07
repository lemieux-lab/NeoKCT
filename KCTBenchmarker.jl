using CairoMakie
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

# Plots individual component sizes of a KCT
function plot_component_sizes(
    n_samples::Int,
    kmer_seq_bytes::Int,
    chunk_ids_bytes::Int,
    n_cids_bytes::Int,
    count_words_bytes::Int,
    bitmap_bytes::Int;
    full_pointer_walkthrough::Bool = false,
    kct = nothing
)
    printstyled("Plotting Component Sizes...\n", color=:green)
    f = Figure()
    xs = ["K-mers", "Indexes", "Δ-Positions", "Counts", "Bitmap"]
    ys = [kmer_seq_bytes, chunk_ids_bytes, n_cids_bytes, count_words_bytes, bitmap_bytes]
    Axis(f[1, 1],
         title    = "$(n_samples) Samples NeoKCT - Component Sizes",
         subtitle = "Total Components Size: $(Base.format_bytes(sum(ys)))" *
                    (full_pointer_walkthrough && kct !== nothing ?
                     " - [Full Pointers Walkthrough: $(Base.format_bytes(Base.summarysize(kct)))]" : ""),
         ylabel   = "Size (Bytes)",
         xticks   = (1:length(xs), xs))
    barplot!(ys, color = ys, strokecolor = :black, strokewidth = 1, bar_labels = ys)
    return f
end

# Plots growth of counts words/bitmap with added samples
function plot_counts_growth(history::Vector{Any})
    printstyled("Plotting counts size growth...\n", color=:green)
    all_samples = [e["samples"] for e in history]
    count_words_hist = [e["count_words_bytes"] for e in history]
    bitmap_hist = [e["bitmap_bytes"] for e in history]
    f = Figure()
    ax = Axis(f[1, 1],
              title = "Count Table Size Over Samples",
              xlabel = "Sample Count",
              ylabel = "Size (Bytes)")
    lines!(ax, all_samples, count_words_hist, label = "Count Words (words)", color = :blue)
    scatter!(ax, all_samples, count_words_hist, color = :blue)
    lines!(ax, all_samples, bitmap_hist, label = "Bitmap (bitmap)", color = :orange)
    scatter!(ax, all_samples, bitmap_hist, color = :orange)
    axislegend(ax, position = :lt)
    return f
end

# Plot growth of seqs/n_cids/cids with added samples
function plot_table_growth(history::Vector{Any})
    printstyled("Plotting table size growth...\n", color=:green)
    all_samples = [e["samples"]         for e in history]
    chunk_ids_hist = [e["chunk_ids_bytes"] for e in history]
    kmer_seq_hist = [e["kmer_seq_bytes"]  for e in history]
    n_cids_hist = [e["n_cids_bytes"]   for e in history]
    f = Figure()
    ax = Axis(f[1, 1],
              title = "CSR Table Sizes Over Samples",
              xlabel = "Sample Count",
              ylabel = "Size (Bytes)")
    lines!(ax, all_samples, chunk_ids_hist, label = "Chunk IDs (flat_cids)", color = :olive)
    scatter!(ax, all_samples, chunk_ids_hist, color = :olive)
    lines!(ax, all_samples, kmer_seq_hist, label = "K-mer Sequences (seqs)", color = :darkgreen)
    scatter!(ax, all_samples, kmer_seq_hist, color = :darkgreen)
    lines!(ax, all_samples, n_cids_hist, label = "Chunk IDs Δ-Positions (n_cids)", color = :lightgreen)
    scatter!(ax, all_samples, n_cids_hist, color = :lightgreen)
    axislegend(ax, position = :lt)
    return f
end

# Plot time between benchmark calls in history (plots individual & total elapsed)
function plot_elapsed_time(history::Vector{Any})
    printstyled("Plotting Elapsed Time...\n", color=:green)
    all_samples = [e["samples"] for e in history]
    timestamps = DateTime.(string.(e["timestamp"]) for e in history)
    delta_samples = all_samples[2:end]
    delta_hours = [Dates.value(timestamps[i] - timestamps[i-1]) / 3_600_000 for i in 2:length(timestamps)]
    total_hours = [Dates.value(t - timestamps[1]) / 3_600_000 for t in timestamps]
    f = Figure(size = (900, 700))
    ax_delta = Axis(f[1, 1],
                    title = "Elapsed Time Between Benchmark Calls",
                    xlabel = "Sample Count",
                    ylabel = "Elapsed Time (hours)")
    barplot!(ax_delta, delta_samples, delta_hours, color = :steelblue)
    ax_total = Axis(f[2, 1],
                    title = "Total Elapsed Time Since First Benchmark Call",
                    subtitle = "Time zero measured at $(all_samples[1]) samples",
                    xlabel = "Sample Count",
                    ylabel = "Total Time (hours)")
    lines!(ax_total, all_samples, total_hours, color = :crimson)
    scatter!(ax_total, all_samples, total_hours, color = :crimson)
    return f
end

# Plot progression of query speed of benchark_size
function plot_query_speed(history::Vector{Any}, benchmark_size::Int)
    printstyled("Plotting Query Speed...\n", color=:green)
    all_samples = [e["samples"]       for e in history]
    query_time_hist = [e["query_time_ms"] for e in history]
    f = Figure()
    ax = Axis(f[1, 1],
              title = "Query Speed Over Samples ($benchmark_size k-mer lookups)",
              xlabel = "Sample Count",
              ylabel = "Query Time (ms)")
    lines!(ax, all_samples, query_time_hist, color = :purple)
    scatter!(ax, all_samples, query_time_hist, color = :purple)
    return f
end

# Plot progression of k-mer unicity
function plot_kmer_cardinality(history::Vector{Any})
    printstyled("Plotting K-mer Cardinality...\n", color=:green)
    all_samples = [e["samples"] for e in history]
    n_kmers_hist = [e["n_kmers"] for e in history]
    f = Figure()
    ax = Axis(f[1, 1],
              title = "Unique K-mer Count Over Samples",
              subtitle = "Saturation indicates diminishing new k-mers per sample",
              xlabel = "Sample Count",
              ylabel = "Unique K-mers")
    lines!(ax, all_samples, n_kmers_hist, color = :teal)
    scatter!(ax, all_samples, n_kmers_hist, color = :teal)
    return f
end

# Plot everything on a KCT
function benchmark_kct(kct::NeoKCT{K, Ab}, benchmark_path::String; full_pointer_walkthrough::Bool=false) where {K, Ab<:Alphabet}
    mkpath(benchmark_path * "sizes_benchmarks/")
    benchmark_file = benchmark_path * "benchmark_data.json"
    history = load_benchmark_history(benchmark_file)

    printstyled("Measuring KCT components...\n", color=:green)
    benchmark_size = 100_000_000
    n_samples = kct.samples.x
    n_kmers = length(kct.seqs)
    kmer_seq_bytes = sizeof(UInt64) * length(kct.seqs)
    chunk_ids_bytes = sizeof(UInt32) * length(kct.flat_cids)
    n_cids_bytes = sizeof(UInt32) * length(kct.n_cids)
    count_words_bytes = Base.summarysize(kct.counts.words)
    bitmap_bytes = Base.summarysize(kct.counts.bitmap)

    k_mers = rand(kct.seqs, benchmark_size)
    t_start = now()
    @showprogress "Benchmarking Query Speed for $benchmark_size k-mers..." for k_mer in k_mers
        findfirst(kct, k_mer)
    end
    query_time_ms = Dates.value(now() - t_start)

    printstyled("Updating Benchmark History...\n", color=:green)
    entry = Dict{String, Any}(
        "samples" => n_samples,
        "timestamp" => string(now()),
        "kmer_seq_bytes" => kmer_seq_bytes,
        "chunk_ids_bytes" => chunk_ids_bytes,
        "n_cids_bytes" => n_cids_bytes,
        "count_words_bytes" => count_words_bytes,
        "bitmap_bytes" => bitmap_bytes,
        "n_kmers" => n_kmers,
        "query_time_ms" => query_time_ms
    )
    save_benchmark_history!(history, entry, benchmark_file)

    f1 = plot_component_sizes(n_samples,
                         kmer_seq_bytes, chunk_ids_bytes, n_cids_bytes, count_words_bytes, bitmap_bytes;
                         full_pointer_walkthrough=full_pointer_walkthrough, kct=kct)

    CairoMakie.save(benchmark_path * "sizes_benchmarks/sizes_benchmark_$(n_samples)_samples.svg", f1)

    length(history) <= 1 && return

    # History dependent plotting
    f2 = plot_counts_growth(history)
    CairoMakie.save(benchmark_path * "counts_benchmark.svg", f2)
    f3 = plot_table_growth(history)
    CairoMakie.save(benchmark_path * "table_benchmark.svg", f3)
    f4 = plot_elapsed_time(history)
    CairoMakie.save(benchmark_path * "elapsed_time_benchmark.svg", f4)
    f5 = plot_query_speed(history, benchmark_size)
    CairoMakie.save(benchmark_path * "query_speed_benchmark.svg", f5)
    f6 = plot_kmer_cardinality(history)
    CairoMakie.save(benchmark_path * "kmer_cardinality_benchmark.svg", f6)

    return
end