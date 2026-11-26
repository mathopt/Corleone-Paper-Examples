# Runs all examples of Section 3, prints a summary of the BenchmarkTools trials, and stores
# them in benchmarks.json (in RESULTS_DIR, default: results/). The MINLP example
# (about 1.5 minutes per solve) only runs with RUN_MINLP=true.
#
# Run:  julia +1.12 --project=. -t 4 scripts/run_all.jl
#       RUN_MINLP=true julia +1.12 --project=. -t 4 scripts/run_all.jl

using BenchmarkTools, Printf
using InteractiveUtils

examples = ["fishing_single_shooting", "fishing_multiple_shooting",
            "fishing_multiple_shooting_sparse_jacobian",
            "fishing_multiple_shooting_sparse_hessian", "fishing_cia", "ude", "oed",
            "oed_cia"]
get(ENV, "RUN_MINLP", "false") == "true" && push!(examples, "fishing_minlp")

results = BenchmarkGroup()
for name in examples
    println("\n=== ", name)
    m = Module(Symbol(:Example_, name))
    Base.include(m, joinpath(@__DIR__, "..", "examples", name * ".jl"))   # each example in its own namespace
    results[name] = Base.invokelatest(getproperty, m, :trial)
end

println()
InteractiveUtils.versioninfo()
@printf("%-28s %8s %12s %12s %12s\n", "example", "samples", "min [s]", "median [s]", "memory [MiB]")
for name in examples
    t = results[name]
    @printf("%-28s %8d %12.4f %12.4f %12.1f\n", name, length(t.times),
            minimum(t).time / 1e9, median(t).time / 1e9, memory(t) / 2^20)
end
BenchmarkTools.save(joinpath(get(ENV, "RESULTS_DIR", mkpath(joinpath(@__DIR__, "..", "results"))), "benchmarks.json"), results)
