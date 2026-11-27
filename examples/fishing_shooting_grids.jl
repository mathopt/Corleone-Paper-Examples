# Section 3.1: influence of the shooting grid and of the initialization of the lifted states on the
# convergence of the Lotka–Volterra fishing problem. Four discretizations of the problem of
# fishing_single_shooting.jl (48 control intervals) are compared:
#   single  single shooting,
#   bad     shooting nodes at t = 3.4, 6.1, 9.75 (off the control grid, as in fishing_multiple_shooting.jl),
#   K11     11 equidistant shooting nodes (every 4th control grid point),
#   K47     47 equidistant shooting nodes (every control grid point).
# All three states are lifted at the nodes. The lifted states are initialized either with the
# initial values of the problem (the default of an annotation without `default`) or with the
# reference x₁ = x₂ = 1 of the objective (and x₃ = 0). All variants use the exact sparse Lagrangian
# Hessian, since the iteration counts of the limited-memory approximation depend on small numerical
# differences in the derivatives (Appendix A). Only iteration counts are reported, no timings.
#
# Run:  julia +1.12 --project=. -t 4 examples/fishing_shooting_grids.jl

using Corleone
using LuxCore, Random
using OrdinaryDiffEqTsit5
using Optimization, OptimizationMOI, Ipopt
using ComponentArrays, ForwardDiff   # ComponentArrays: default parameter vectorizer
import DifferentiationInterface, SparseMatrixColorings   # enable sparse derivatives
using Printf

function lotka!(dx, x, p, t)
    dx[1] = x[1] - p[1] * x[1] * x[2] - 0.4 * p[2] * x[1]
    dx[2] = -x[2] + p[3] * x[1] * x[2] - 0.2 * p[2] * x[2]
    dx[3] = (x[1] - 1.0)^2 + (x[2] - 1.0)^2
    return nothing
end
prob = ODEProblem(lotka!, [0.5, 0.7, 0.0], (0.0, 12.0), [1.0, 0.0, 1.0];
                  abstol = 1e-9, reltol = 1e-9)
grid = collect(0.25:0.25:11.75)              # 48 control intervals
reference = [1.0, 1.0, 0.0]                  # x₁ = x₂ = 1 of the objective, no accumulated cost

grids = (single = Float64[], bad = [3.4, 6.1, 9.75],
         K11 = collect(1.0:1.0:11.0), K47 = collect(0.25:0.25:11.75))

function fishing_problem(nodes; init = :initial)
    control = Tunable(2; kind = Val(:parameter), default = 1.0,
                      lowerbound = 0.0, upperbound = 1.0) => grid
    fixed = Constant([3, 1]; kind = Val(:parameter)) => nothing
    lifted = (init == :reference ?
              Tunable(i; kind = Val(:variable), default = reference[i],
                      lowerbound = 0.0, upperbound = 2.0) :
              Tunable(i; kind = Val(:variable), lowerbound = 0.0, upperbound = 2.0)
              for i in 1:3)
    layer = isempty(nodes) ?
            ShootingLayer(prob, Tsit5(), control, fixed; controls = [2], quadratures = [3]) :
            ShootingLayer(prob, Tsit5(), control, fixed, (l => nodes for l in lifted)...;
                          ensemblealg = EnsembleThreads(), controls = [2], quadratures = [3])
    ps, st = LuxCore.setup(MersenneTwister(42), layer)
    ps, st = add_parameter_segments!(layer, ps, st)
    dynprob = CorleoneDynamicOptProblem(layer, :(u₃(12.0)); ps, st)
    opt = OptimizationProblem(dynprob, AutoForwardDiff(); sparse_constraints = !isempty(nodes),
                              sparse_hessian = true)
    return (; layer, ps, st, dynprob, opt)
end

solve_exact(opt) = solve(opt, Ipopt.Optimizer(); max_iter = 1000, print_level = 0, tol = 1e-6)

results = []
for name in keys(grids), init in (name == :single ? (:none,) : (:initial, :reference))
    case = fishing_problem(grids[name]; init)
    sol = solve_exact(case.opt)
    @assert successful_retcode(sol)
    push!(results, (; name, init, case..., sol, vars = length(sol.u), iters = sol.stats.iterations))
    @printf("%-7s %-10s variables = %3d  iterations = %3d  objective = %.6f\n",
            name, init, length(sol.u), sol.stats.iterations, sol.objective)
end
println("threads = ", Threads.nthreads())
