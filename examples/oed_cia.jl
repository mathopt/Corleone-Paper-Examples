# Section 3.3: binary design from the relaxed OED solution by sum-up rounding (combinatorial
# integral approximation). Solves the relaxed problem of oed.jl, rounds the measurement weights
# w1, w2 on their 48 intervals, keeps the relaxed fishing control, and evaluates the A-criterion
# and the measurement budgets of the rounded design.
#
# The benchmark times the rounding and the evaluation of the rounded design; the relaxed solve is
# timed by oed.jl, which this script includes.
#
# Run:  julia +1.12 --project=. examples/oed_cia.jl

Base.include(@__MODULE__, joinpath(@__DIR__, "oed.jl"))   # defines oed, ps, st, opt, sol, reconstruct, traj, ...

using ComponentArrays: labels

# Sum-up rounding on an equidistant grid: w_bin[k] = 1 if the accumulated relaxed weight
# exceeds the accumulated binary weight by at least half an interval.
function sum_up_rounding(w)
    wbin = zeros(length(w))
    acc = 0.0
    for k in eachindex(w)
        acc += w[k]
        wbin[k] = acc - sum(wbin[1:k-1]) >= 0.5 ? 1.0 : 0.0
    end
    return wbin
end

# Positions of the measurement weights in the flat parameter vector; their ComponentArray labels
# read parameter.parameters.p₁.parameter.w1.ps[k].
psol = reconstruct(sol.u)
weight_indices(p, name) = findall(l -> occursin(".$(name).", l), labels(p))

for name in (:w1, :w2)
    w = sol.u[weight_indices(psol, name)]
    println(name, ": fractional relaxed weights (1e-3 < w < 1 - 1e-3) = ",
            count(x -> 1e-3 < x < 1 - 1e-3, w), " of ", length(w),
            ", max. distance to 0 or 1 = ", maximum(abs.(w .- round.(w))))
end

function round_design(u)
    u_rounded = copy(u)
    for name in (:w1, :w2)
        idx = weight_indices(psol, name)
        u_rounded[idx] .= sum_up_rounding(u[idx])
    end
    traj_rounded, _ = oed(nothing, reconstruct(u_rounded), st)
    return traj_rounded, ACriterion()(fisher_information(oed, traj_rounded))
end

traj_rounded, criterion_rounded = round_design(sol.u)
println("Relaxed:  A-criterion = ", sol.objective, ", budgets = ", sampling_sums(oed, traj))
println("Rounded:  A-criterion = ", criterion_rounded,
        ", budgets = ", sampling_sums(oed, traj_rounded))
println("Note: sum-up rounding can exceed the budget of 4 by at most half an interval (0.125).")

trial = @benchmark round_design($(sol.u)) samples = 10 evals = 1 seconds = 60
display(trial)
println()
