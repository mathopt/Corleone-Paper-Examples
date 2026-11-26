# Section 3.1: Lotka–Volterra fishing problem with multiple shooting. All three states are
# lifted at t = 3.4, 6.1, 9.75; the shooting blocks are solved with a threaded ensemble.
#
# Run:  julia +1.12 --project=. -t 4 examples/fishing_multiple_shooting.jl

using Corleone
using LuxCore, Random
using OrdinaryDiffEqTsit5
using Optimization, OptimizationMOI, Ipopt
using ComponentArrays, ForwardDiff   # ComponentArrays: default parameter vectorizer
using BenchmarkTools

function lotka!(dx, x, p, t)
    dx[1] = x[1] - p[1] * x[1] * x[2] - 0.4 * p[2] * x[1]
    dx[2] = -x[2] + p[3] * x[1] * x[2] - 0.2 * p[2] * x[2]
    dx[3] = (x[1] - 1.0)^2 + (x[2] - 1.0)^2
    return nothing
end
prob = ODEProblem(lotka!, [0.5, 0.7, 0.0], (0.0, 12.0), [1.0, 0.0, 1.0];
                  abstol = 1e-9, reltol = 1e-9)
grid = collect(0.25:0.25:11.75)              # 48 control intervals
control = Tunable(2; kind = Val(:parameter), default = 1.0,
                  lowerbound = 0.0, upperbound = 1.0) => grid
fixed = Constant([3, 1]; kind = Val(:parameter)) => nothing
shooting_points = [3.4, 6.1, 9.75]
lifted = (Tunable(i; kind = Val(:variable), lowerbound = 0.0, upperbound = 2.0) => shooting_points
          for i in 1:3)
layer = ShootingLayer(prob, Tsit5(), control, fixed, lifted...;
                      ensemblealg = EnsembleThreads(), controls = [2], quadratures = [3])
ps, st = LuxCore.setup(MersenneTwister(42), layer)
ps, st = add_parameter_segments!(layer, ps, st)
dynprob = CorleoneDynamicOptProblem(layer, :(u₃(12.0)); ps, st)
opt = OptimizationProblem(dynprob, AutoForwardDiff())

run_example() = solve(opt, Ipopt.Optimizer(); max_iter = 1000,
                      hessian_approximation = "limited-memory", print_level = 0, tol = 1e-6)

sol = run_example()
@assert successful_retcode(sol)
_, reconstruct = pack(default_vectorizer(), dynprob.ps)
traj, _ = layer(nothing, reconstruct(sol.u), opt.p)
println("Multiple shooting: objective = ", sol.objective,
        ", max. continuity residual = ", maximum(abs, shooting_constraints(traj)),
        ", iterations = ", sol.stats.iterations, ", variables = ", length(sol.u),
        ", threads = ", Threads.nthreads())

trial = @benchmark run_example() samples = 10 evals = 1 seconds = 60
display(trial)
println()
