# Figure of Section 3.2 (figures/ude-fit.pdf), full page width. Left: trained universal
# differential equation, noisy observations, and the true model. Right: learned interaction
# N_θ(x)² along the fitted trajectory and the true interaction x₁x₂ along the true trajectory.
# Runs examples/ude.jl (fit and timing) and plots the result.
#
# Run:  julia +1.12 --project=. figures/ude_figure.jl   (output: FIGURE_DIR, default: results/)

include(joinpath(@__DIR__, "..", "examples", "ude.jl"))   # defines prob, weights, model, st, sol, data, truth, ...

using CairoMakie, LaTeXStrings
using SymbolicIndexingInterface: parameter_values_at_time

# Fitted parameters, evaluated on a fine time grid with the same annotations as the fitted layer.
tfine = 0.0:0.01:5.0
_, reconstruct = pack(default_vectorizer(), model.ps)
fine = ShootingLayer(remake(prob; saveat = tfine), Tsit5(), weights)
fitted, _ = fine(nothing, reconstruct(sol.u), st)
X = reduce(hcat, fitted.u)
θ = parameter_values_at_time(fitted, fitted, 0.0)   # trained network weights
learned = [only(first(Lux.apply(network, x, rebuild(θ), nn_st)))^2 for x in fitted.u]

truth_fine = solve(remake(truth.prob; saveat = tfine), Tsit5(); abstol = 1e-9, reltol = 1e-9)
T = reduce(hcat, truth_fine.u)

colors = Makie.wong_colors()
theme = merge(theme_latexfonts(), Theme(fontsize = 22))
with_theme(theme) do
    fig = Figure(; size = (1200, 420))
    states = Makie.Axis(fig[1, 1]; xlabel = L"t", ylabel = L"x(t)")
    for i in 1:2
        lines!(states, fitted.t, X[i, :]; color = colors[i], linewidth = 2)
        scatter!(states, times, data[i, :]; color = colors[i], markersize = 8)
        lines!(states, truth_fine.t, T[i, :]; color = :black, linestyle = :dash)
    end
    interaction = Makie.Axis(fig[1, 2]; xlabel = L"t", ylabel = "Interaction")
    lines!(interaction, fitted.t, learned; color = colors[3], linewidth = 2)
    lines!(interaction, truth_fine.t, T[1, :] .* T[2, :]; color = :black, linestyle = :dash)

    elements = [
        LineElement(; color = colors[1], linewidth = 2), LineElement(; color = colors[2], linewidth = 2),
        MarkerElement(; color = colors[1], marker = :circle, markersize = 8),
        MarkerElement(; color = colors[2], marker = :circle, markersize = 8),
        LineElement(; color = colors[3], linewidth = 2),
        LineElement(; color = :black, linestyle = :dash, linewidth = 2),
    ]
    labels = [L"x_1", L"x_2", "Prey observations", "Predator observations",
              L"\text{Learned } N_\theta(x)^2", L"\text{Groundtruth model, } x_1 x_2"]
    Legend(fig[2, 1:2], elements, labels; orientation = :horizontal, nbanks = 2,
           framevisible = false, colgap = 20, rowgap = 2)
    out = joinpath(get(ENV, "FIGURE_DIR", mkpath(joinpath(@__DIR__, "..", "results"))), "ude-fit.pdf")
    save(out, fig)
    save(replace(out, ".pdf" => ".png"), fig)
    println("saved ", normpath(out))
end
