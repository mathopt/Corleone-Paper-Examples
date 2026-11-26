# Section 3.1, "Integer controls": Listing 1 with u(t) ∈ {0, 1}, solved with the
# branch-and-bound solver Juniper and Ipopt. Integrality is declared on the control
# annotation (`integer = true`); OptimizationProblem derives the `int` mask from it.
# A single solve takes about 1.5 minutes; 10 samples are taken as for the other
# examples (set BENCHMARK_SAMPLES to change the number).
#
# Run:  julia +1.12 --project=. examples/fishing_minlp.jl

using Corleone
using OrdinaryDiffEqTsit5
using Optimization, OptimizationMOI, Ipopt, Juniper
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
                  lowerbound = 0.0, upperbound = 1.0, integer = true) => grid
fixed = Constant([3, 1]; kind = Val(:parameter)) => nothing
layer = ShootingLayer(prob, Tsit5(), control, fixed;
                      controls = [2], quadratures = [3])
dynprob = CorleoneDynamicOptProblem(layer, :(u₃(12.0)))
opt = OptimizationProblem(dynprob, AutoForwardDiff())
@assert opt.int == trues(48)                 # mask derived from the annotation

ipopt = OptimizationMOI.MOI.OptimizerWithAttributes(Ipopt.Optimizer, "print_level" => 0)
run_example() = solve(opt, Juniper.Optimizer(); nl_solver = ipopt,
                      log_levels = Symbol[])

sol = run_example()
@assert all(u -> isapprox(u, 0; atol = 1e-6) || isapprox(u, 1; atol = 1e-6), sol.u)
println("MINLP: objective = ", sol.objective, ", retcode = ", sol.retcode,
        ", branches = ", sol.original.inner.nbranches, ", u = ", round.(Int, sol.u))

samples = parse(Int, get(ENV, "BENCHMARK_SAMPLES", "10"))
trial = @benchmark run_example() samples = samples evals = 1 seconds = 3600
display(trial)
println()
