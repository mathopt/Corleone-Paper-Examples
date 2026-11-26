# Section 3.1: integer fishing control from the relaxed solution by sum-up rounding (combinatorial
# integral approximation). Solves the relaxed single-shooting problem of
# fishing_single_shooting.jl, rounds the 48 control values to {0, 1}, and evaluates the objective
# of the rounded control. Compare with the MINLP solution of fishing_minlp.jl.
#
# The benchmark times the rounding and the evaluation of the rounded control; the relaxed solve is
# timed by fishing_single_shooting.jl, which this script includes.
#
# Run:  julia +1.12 --project=. examples/fishing_cia.jl

Base.include(@__MODULE__, joinpath(@__DIR__, "fishing_single_shooting.jl"))   # defines opt, sol, ...

# Sum-up rounding on an equidistant grid: w_bin[k] = 1 if the accumulated relaxed control
# exceeds the accumulated binary control by at least half an interval.
function sum_up_rounding(w)
    wbin = zeros(length(w))
    acc = 0.0
    for k in eachindex(w)
        acc += w[k]
        wbin[k] = acc - sum(wbin[1:k-1]) >= 0.5 ? 1.0 : 0.0
    end
    return wbin
end

# Maximum deviation of the accumulated controls (integral approximation error), Δt = 0.25.
integral_deviation(w, wbin) = 0.25 * maximum(abs, cumsum(w .- wbin))

objective(u) = opt.f(u, opt.p)                # x₃(12) of the control values u

function round_control(u)
    ubin = sum_up_rounding(u)
    return ubin, objective(ubin)
end

u = collect(sol.u)                            # all 48 optimization variables are control values
@assert objective(u) ≈ sol.objective
ubin, objective_rounded = round_control(u)
println("Relaxed: objective = ", sol.objective, ", fractional control values (1e-3 < u < 1 - 1e-3) = ",
        count(x -> 1e-3 < x < 1 - 1e-3, u), " of ", length(u))
println("Rounded: objective = ", objective_rounded, ", switches = ", count(!iszero, diff(ubin)),
        ", integral deviation = ", integral_deviation(u, ubin),
        ", gap to relaxed = ", objective_rounded - sol.objective)

trial = @benchmark round_control($u) samples = 10 evals = 1 seconds = 60
display(trial)
println()
