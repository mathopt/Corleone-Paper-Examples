# Section 3.2, Listing 2: universal differential equation for the Lotka–Volterra model,
# trained on 51 noisy observations of both populations with CorleoneInverse
# (weighted least squares, Levenberg–Marquardt from NonlinearSolve).
#
# Run:  julia +1.12 --project=. examples/ude.jl

using CorleoneInverse                        # re-exports Corleone
using Lux, Random, ComponentArrays
using OrdinaryDiffEqTsit5
using NonlinearSolve
using ForwardDiff
using BenchmarkTools

# Synthetic data from the true model with interaction x₁x₂ and noise σ = 0.05.
rng = MersenneTwister(42)
times = collect(0.0:0.1:5.0)
σ = 0.05
true_lotka(x, p, t) = [x[1] - x[1] * x[2], x[1] * x[2] - x[2]]
truth = solve(ODEProblem(true_lotka, [0.8, 1.2], (0.0, 5.0)), Tsit5();
              saveat = times, abstol = 1e-9, reltol = 1e-9)
data = Array(truth) .+ σ .* randn(rng, 2, length(times))

# Network with one hidden layer of five tanh units (21 weights), small initial weights.
network = Lux.Chain(Lux.Dense(2 => 5, tanh), Lux.Dense(5 => 1))
nn_ps, nn_st = Lux.setup(rng, network)
nn_ps = ComponentArray{Float64}(nn_ps)
nn_ps .*= 0.1
nn_ps.layer_2.bias .= 0.7
p0, rebuild = pack(default_vectorizer(), nn_ps)
nn_st = Lux.testmode(nn_st)
function hybrid(x, p, t)
    out, _ = Lux.apply(network, x, rebuild(p), nn_st)
    interaction = only(out)^2
    return [x[1] - interaction, interaction - x[2]]
end
prob = ODEProblem(hybrid, [0.8, 1.2], (0.0, 5.0), p0;
                  abstol = 1e-8, reltol = 1e-8)
weights = Tunable(collect(eachindex(p0)); kind = Val(:parameter),
                  lowerbound = -3.0, upperbound = 3.0) => nothing
layer = ShootingLayer(prob, Tsit5(), weights)
ps, st = Lux.setup(rng, layer)
experiment = Experiment(1 => (times, data[1, :], σ),   # observations of state 1 (prey)
                        2 => (times, data[2, :], σ))   # and state 2 (predators)
model = CorleoneDynamicOptProblem(layer, experiment; ps, st)
nlls = NonlinearLeastSquaresProblem(model)

run_example() = solve(nlls, LevenbergMarquardt(); maxiters = 100)

# Residuals are weighted by 1/σ, so σ² times their mean square is the mean squared error.
mse(residuals) = σ^2 * sum(abs2, residuals) / length(residuals)
sol = run_example()
@assert all(isfinite, sol.resid)
println("UDE: weights = ", length(p0), ", observations = ", length(sol.resid),
        ", initial MSE = ", σ^2 * mean_squared_residual(model),
        ", final MSE = ", mse(sol.resid), ", iterations = ", sol.stats.nsteps,
        ", retcode = ", sol.retcode)

trial = @benchmark run_example() samples = 10 evals = 1 seconds = 60
display(trial)
println()
