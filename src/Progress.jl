## Progress ##

# Shared progress reporter for the paper's result-generating scripts (see PAPER_TODO.md's
# conventions section): anything over ~1e7 items or ~1 minute logs count, rate, ETA and RSS
# through this instead of @showprogress, since @showprogress's own overhead shows up in
# timing loops that use time_ns() around the whole call.

# println + flush: stdout is block-buffered when redirected to a file, so without the flush
# a multi-hour run shows nothing in the log until the buffer fills or the process exits.
# Locked so concurrent callers (StreamBuild's shard-parallel loops, _Prog's tick!) don't
# interleave half-lines.
const _LOG_LOCK = ReentrantLock()
_log(args...) = lock(_LOG_LOCK) do
    println("[", Dates.format(now(), "HH:MM:SS"), "] ", args...); flush(stdout)
end

# Current resident set size in GB (Linux). Sys.maxrss() only gives the peak, not the
# current value, which is what a "is this run about to OOM" check needs.
_rss_gb() = parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 / 1e9

mutable struct _Prog
    total::Int
    done::Threads.Atomic{Int}
    t0::UInt64
    last::Threads.Atomic{UInt64}
    every::Float64
    label::String
end

_Prog(total::Integer, label::String; every::Real=30) =
    _Prog(Int(total), Threads.Atomic{Int}(0), time_ns(), Threads.Atomic{UInt64}(time_ns()),
          Float64(every), label)

# Call from every worker after finishing n items (default 1). Thread-safe: only the thread
# that wins the atomic swap on `last` prints, so concurrent callers never interleave a line.
function tick!(p::_Prog, n::Integer=1)
    d = Threads.atomic_add!(p.done, Int(n)) + Int(n)
    t = time_ns()
    (t - p.last[]) / 1e9 < p.every && d < p.total && return
    Threads.atomic_xchg!(p.last, t)
    el = (t - p.t0) / 1e9
    rate = d / el
    eta = (p.total - d) / max(rate, 1e-9)
    _log(p.label, ": ", d, "/", p.total, " (", round(100d / p.total; digits=2), "%), ",
         round(rate / 1e6; digits=3), " M/s, ETA ", round(eta / 3600; digits=2), " h, RSS ",
         round(_rss_gb(); digits=1), " GB")
end
