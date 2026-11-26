# Section 3.3, Listing 3: optimal experimental design for the interaction coefficients of
# the Lotka–Volterra fishing model with CorleoneOED (A-criterion, budget 4 per population).
#
# Run:  julia +1.12 --project=. examples/oed.jl

using Corleone, CorleoneOED
using LuxCore, Random
using OrdinaryDiffEqTsit5
using Optimization, OptimizationMOI, Ipopt
using ComponentArrays, ForwardDiff   # ComponentArrays: default parameter vectorizer
using BenchmarkTools

# p[1]: fishing intensity, p[2], p[3]: estimated interaction coefficients.
function lotka!(dx, x, p, t)
    dx[1] = x[1] - p[2] * x[1] * x[2] - 0.4 * p[1] * x[1]
    dx[2] = -x[2] + p[3] * x[1] * x[2] - 0.2 * p[1] * x[2]
    return nothing
end
prob = ODEProblem(lotka!, [0.5, 0.7], (0.0, 12.0), [0.0, 1.0, 1.0])
fishing = Tunable(1; default = 0.0, lowerbound = 0.0,
                  upperbound = 1.0) => collect(0.25:0.25:11.75)
grid = collect(0.0:0.25:11.75)
measurements = (
    ContinuousMeasurement(:w1, grid, (x, p, t) -> x[1]),
    ContinuousMeasurement(:w2, grid, (x, p, t) -> x[2]),
)
oed = OEDLayer(prob, Tsit5(), fishing; params = [2, 3], measurements,
               controls = [1])
ps, st = LuxCore.setup(Xoshiro(1), oed)
model = CorleoneDynamicOptProblem(oed, ACriterion(); ps, st,
                                  M = [4.0, 4.0])
opt = OptimizationProblem(model, AutoForwardDiff())

# Baseline: no fishing, both populations observed uniformly with the same budget of 4.
# The initial design observes both with weight 1 on [0, 12] (budget 12); the information
# is linear in the weights, so the uniform budget-4 design has information F / 3.
traj0, _ = oed(nothing, ps, st)
@assert all(isapprox.(sampling_sums(oed, traj0), 12.0; rtol = 1e-8))
baseline = ACriterion()(fisher_information(oed, traj0) ./ 3)

run_example() = solve(opt, Ipopt.Optimizer(); max_iter = 100,
                      hessian_approximation = "limited-memory", print_level = 0, tol = 1e-6)

sol = run_example()
@assert successful_retcode(sol)
_, reconstruct = pack(default_vectorizer(), ps)
traj, _ = oed(nothing, reconstruct(sol.u), st)
println("OED: A-criterion = ", sol.objective,
        ", budgets = ", sampling_sums(oed, traj), ", iterations = ", sol.stats.iterations,
        ", baseline = ", baseline, ", improvement = ", baseline / sol.objective)

trial = @benchmark run_example() samples = 10 evals = 1 seconds = 60
display(trial)
println()
