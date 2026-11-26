# Section 3.1, Listing 1: Lotka–Volterra fishing problem with single shooting
# (48 control intervals, Ipopt with limited-memory Hessian approximation).
#
# Run:  julia +1.12 --project=. examples/fishing_single_shooting.jl

using Corleone
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
layer = ShootingLayer(prob, Tsit5(), control, fixed;
                      controls = [2], quadratures = [3])
dynprob = CorleoneDynamicOptProblem(layer, :(u₃(12.0)))
opt = OptimizationProblem(dynprob, AutoForwardDiff())

run_example() = solve(opt, Ipopt.Optimizer();
                      hessian_approximation = "limited-memory", print_level = 0, tol = 1e-6)

sol = run_example()
@assert successful_retcode(sol)
println("Single shooting: objective = ", sol.objective,
        ", iterations = ", sol.stats.iterations, ", variables = ", length(sol.u))

trial = @benchmark run_example() samples = 10 evals = 1 seconds = 60
display(trial)
println()
