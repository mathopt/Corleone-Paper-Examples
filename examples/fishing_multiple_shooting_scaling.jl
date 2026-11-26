# Section 3.1: scaling of multiple shooting for the Lotka–Volterra fishing problem. Solves the
# problem of fishing_multiple_shooting.jl with K = 0 (single shooting), 1, 3, 5, 11, 23, and 47
# equidistant shooting nodes (all three states lifted, 48 control intervals as before), with a
# dense and a sparse constraint Jacobian, and reports iterations and the median solve time.
# The shooting segments run in a threaded ensemble; run the script with different numbers of
# threads to see the effect of parallelization (scaling.sh loops over 1, 2, 4, and 8 threads).
#
# Run:  julia +1.12 --project=. -t 4 examples/fishing_multiple_shooting_scaling.jl
#       SCALING_NODES=0,3,11 SCALING_SAMPLES=3 julia +1.12 --project=. -t 4 examples/fishing_multiple_shooting_scaling.jl

using Corleone
using LuxCore, Random
using OrdinaryDiffEqTsit5
using Optimization, OptimizationMOI, Ipopt
using ComponentArrays, ForwardDiff   # ComponentArrays: default parameter vectorizer
import DifferentiationInterface, SparseMatrixColorings   # enable sparse derivatives
using BenchmarkTools, Printf

function lotka!(dx, x, p, t)
    dx[1] = x[1] - p[1] * x[1] * x[2] - 0.4 * p[2] * x[1]
    dx[2] = -x[2] + p[3] * x[1] * x[2] - 0.2 * p[2] * x[2]
    dx[3] = (x[1] - 1.0)^2 + (x[2] - 1.0)^2
    return nothing
end
prob = ODEProblem(lotka!, [0.5, 0.7, 0.0], (0.0, 12.0), [1.0, 0.0, 1.0];
                  abstol = 1e-9, reltol = 1e-9)
grid = collect(0.25:0.25:11.75)              # 48 control intervals

function fishing_problem(K; sparse)
    control = Tunable(2; kind = Val(:parameter), default = 1.0,
                      lowerbound = 0.0, upperbound = 1.0) => grid
    fixed = Constant([3, 1]; kind = Val(:parameter)) => nothing
    nodes = collect(range(0.0, 12.0; length = K + 2)[2:(end - 1)])   # K + 1 equal segments
    lifted = (Tunable(i; kind = Val(:variable), lowerbound = 0.0, upperbound = 2.0) => nodes
              for i in 1:3 if K > 0)
    layer = ShootingLayer(prob, Tsit5(), control, fixed, lifted...;
                          ensemblealg = EnsembleThreads(), controls = [2], quadratures = [3])
    ps, st = LuxCore.setup(MersenneTwister(42), layer)
    ps, st = add_parameter_segments!(layer, ps, st)
    dynprob = CorleoneDynamicOptProblem(layer, :(u₃(12.0)); ps, st)
    return OptimizationProblem(dynprob, AutoForwardDiff(); sparse_constraints = sparse && K > 0)
end

solve_problem(opt) = solve(opt, Ipopt.Optimizer(); max_iter = 1000,
                           hessian_approximation = "limited-memory", print_level = 0, tol = 1e-6)

nodes_list = parse.(Int, split(get(ENV, "SCALING_NODES", "0,1,3,5,11,23,47"), ","))
samples = parse(Int, get(ENV, "SCALING_SAMPLES", "5"))
rows = []
for K in nodes_list, sparse in (K == 0 ? (false,) : (false, true))
    opt = fishing_problem(K; sparse)
    sol = solve_problem(opt)                 # compilation and check
    trial = @benchmark solve_problem($opt) samples = samples evals = 1 seconds = 300
    push!(rows, (; K, sparse, vars = length(sol.u), iters = sol.stats.iterations,
                 objective = sol.objective, ok = successful_retcode(sol),
                 time = median(trial).time / 1e9))
    r = rows[end]
    @printf("K = %2d  %-6s  variables = %3d  iterations = %4d  objective = %.6f  %s  median = %.3f s\n",
            K, sparse ? "sparse" : "dense", r.vars, r.iters, r.objective, r.ok ? "ok" : "FAILED", r.time)
end

println("\nthreads = ", Threads.nthreads())
@printf("%5s %8s %10s %8s %12s %12s\n", "nodes", "jacobian", "variables", "iter.", "objective", "median [s]")
for r in rows
    @printf("%5d %8s %10d %8d %12.6f %12.3f\n", r.K, r.K == 0 ? "--" : (r.sparse ? "sparse" : "dense"),
            r.vars, r.iters, r.objective, r.time)
end
