using CairoMakie
using JSON

# TODO: Break down plotting into multiple functions (one per plot + reading/updating the benchmark_file)
function benchmark_kct(kct::NeoKCT{K, Ab}, benchmark_path::String; full_pointer_walkthrough::Bool=false) where {K, Ab<:Alphabet}
    mkpath(benchmark_path * "sizes_benchmarks/")
    benchmark_file = benchmark_path * "benchmark_data.json"
    history = isfile(benchmark_file) ? JSON.parsefile(benchmark_file) : []

    printstyled("Measuring KCT components...\n", color=:green)
    n_samples = kct.samples.x
    n_kmers = length(kct.seqs)
    kmer_seq_bytes = sizeof(UInt64) * length(kct.seqs)
    chunk_ids_bytes = sizeof(UInt32) * length(kct.flat_cids)
    offsets_bytes = sizeof(UInt32) * length(kct.offsets)
    count_words_bytes = Base.summarysize(kct.counts.words)
    bitmap_bytes = Base.summarysize(kct.counts.bitmap)

    # Query speed benchmark
    benchmark_size = 100_000_000
    k_mers = rand(kct.seqs, benchmark_size)
    t_start = now()
    @showprogress "Benchmarking Query Speed for $benchmark_size k-mers..." for k_mer in k_mers
        findfirst(kct, k_mer)
    end
    query_time_ms = Dates.value(now() - t_start)

    printstyled("Updating Benchark History...\n", color=:green)
    entry = Dict{String, Any}(
        "samples" => n_samples,
        "timestamp" => string(now()),
        "kmer_seq_bytes" => kmer_seq_bytes,
        "chunk_ids_bytes" => chunk_ids_bytes,
        "offsets_bytes" => offsets_bytes,
        "count_words_bytes" => count_words_bytes,
        "bitmap_bytes" => bitmap_bytes,
        "n_kmers" => n_kmers,
        "query_time_ms" => query_time_ms
    )
    push!(history, entry)

    mkpath(benchmark_path)
    open(benchmark_file, "w") do io
        JSON.print(io, history, 2)
    end

    # --- Plot 1: Component sizes bar chart ---
    printstyled("Plotting Component Sizes...\n", color=:green)
    f = Figure()
    xs = ["K-mers", "Indexes", "Counts", "Bitmap"]
    ys = [kmer_seq_bytes, chunk_ids_bytes, count_words_bytes, bitmap_bytes]
    Axis(f[1, 1],
         title = "$(n_samples) Samples NeoKCT - Component Sizes",
         subtitle = "Total Components Size: $(Base.format_bytes(sum(ys)))" * (full_pointer_walkthrough ? " - [Full Pointers Walkthrough: $(Base.format_bytes(Base.summarysize(kct)))]" : ""),
         ylabel = "Size (Bytes)",
         xticks = (1:length(xs), xs))
    barplot!(ys, color = ys, strokecolor = :black, strokewidth = 1, bar_labels = ys)
    CairoMakie.save(benchmark_path * "sizes_benchmarks/sizes_benchmark_$(n_samples)_samples.svg", f)

    length(history) <= 1 && return

    # Extract history arrays
    all_samples = [e["samples"] for e in history]
    count_words_hist = [e["count_words_bytes"] for e in history]
    bitmap_hist = [e["bitmap_bytes"] for e in history]
    chunk_ids_hist = [e["chunk_ids_bytes"] for e in history]
    kmer_seq_hist = [e["kmer_seq_bytes"] for e in history]
    offsets_hist = [e["offsets_bytes"] for e in history]
    n_kmers_hist = [e["n_kmers"] for e in history]
    query_time_hist = [e["query_time_ms"] for e in history]
    timestamps = DateTime.(string.(e["timestamp"]) for e in history)

    # --- Plot 2: Count table size over samples ---
    printstyled("Plotting counts size growth...\n", color=:green)
    f2 = Figure()
    ax2 = Axis(f2[1, 1],
               title  = "Count Table Size Over Samples",
               xlabel = "Sample Count",
               ylabel = "Size (Bytes)")
    lines!(ax2, all_samples, count_words_hist, label = "Count Words (words)", color = :blue)
    scatter!(ax2, all_samples, count_words_hist, color = :blue)
    lines!(ax2, all_samples, bitmap_hist, label = "Bitmap (bitmap)", color = :orange)
    scatter!(ax2, all_samples, bitmap_hist, color = :orange)
    axislegend(ax2, position = :lt)
    CairoMakie.save(benchmark_path * "counts_benchmark.svg", f2)

    # --- Plot 3: Index & sequence pointer sizes over samples ---
    printstyled("Plotting table size growth...\n", color=:green)
    f3 = Figure()
    ax3 = Axis(f3[1, 1],
               title  = "CSR Table Sizes Over Samples",
               xlabel = "Sample Count",
               ylabel = "Size (Bytes)")
    lines!(ax3, all_samples, chunk_ids_hist, label = "Chunk IDs (flat_cids)", color = :green)
    scatter!(ax3, all_samples, chunk_ids_hist, color = :green)
    lines!(ax3, all_samples, kmer_seq_hist, label = "K-mer Sequences (seqs)", color = :darkgreen)
    scatter!(ax3, all_samples, kmer_seq_hist, color = :darkgreen)
    lines!(ax3, all_samples, offsets_hist, label = "CSR offsets (offsets_bytes)", color = :lightgreen)
    scatter!(ax3, all_samples, offsets_hist, color = :lightgreen)
    axislegend(ax3, position = :lt)
    CairoMakie.save(benchmark_path * "table_benchmark.svg", f3)

    # --- Plot 4: Elapsed time (delta + total) ---
    printstyled("Plotting Elapsed Time...\n", color=:green)
    delta_samples = all_samples[2:end]
    delta_hours = [Dates.value(timestamps[i] - timestamps[i-1]) / 3_600_000 for i in 2:length(timestamps)]
    total_hours = [Dates.value(t - timestamps[1]) / 3_600_000 for t in timestamps]

    f4 = Figure(size = (900, 700))
    ax4a = Axis(f4[1, 1],
                title = "Elapsed Time Between Benchmark Calls",
                xlabel = "Sample Count",
                ylabel = "Elapsed Time (hours)")
    barplot!(ax4a, delta_samples, delta_hours, color = :steelblue)

    ax4b = Axis(f4[2, 1],
                title = "Total Elapsed Time Since First Benchmark Call",
                subtitle = "Time zero measured at $(all_samples[1]) samples",
                xlabel = "Sample Count",
                ylabel = "Total Time (hours)")
    lines!(ax4b, all_samples, total_hours, color = :crimson)
    scatter!(ax4b, all_samples, total_hours, color = :crimson)
    CairoMakie.save(benchmark_path * "elapsed_time_benchmark.svg", f4)

    # --- Plot 5: Query speed over samples ---
    printstyled("Plotting Query Speed...\n", color=:green)
    f5 = Figure()
    ax5 = Axis(f5[1, 1],
               title = "Query Speed Over Samples ($benchmark_size k-mer lookups)",
               xlabel = "Sample Count",
               ylabel = "Query Time (ms)")
    lines!(ax5, all_samples, query_time_hist, color = :purple)
    scatter!(ax5, all_samples, query_time_hist, color = :purple)
    CairoMakie.save(benchmark_path * "query_speed_benchmark.svg", f5)

    # --- Plot 6 K-mer cardinality growth ---
    printstyled("Plotting K-mer Cardinality...\n", color=:green)
    f6 = Figure()
    ax6 = Axis(f6[1, 1],
               title = "Unique K-mer Count Over Samples",
               subtitle = "Saturation indicates diminishing new k-mers per sample",
               xlabel = "Sample Count",
               ylabel = "Unique K-mers")
    lines!(ax6, all_samples, n_kmers_hist, color = :teal)
    scatter!(ax6, all_samples, n_kmers_hist, color = :teal)
    CairoMakie.save(benchmark_path * "kmer_cardinality_benchmark.svg", f6)

    return
end