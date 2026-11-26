# Figure of Section 3.3 (figures/oed-optimal-design.pdf), full page width. Left: populations
# x₁, x₂ under the optimized design. Right: fishing intensity u and measurement weights w₁, w₂.
# Runs examples/oed.jl (design and timing) and plots the result.
#
# Run:  julia +1.12 --project=. figures/oed_figure.jl   (output: FIGURE_DIR, default: results/)

include(joinpath(@__DIR__, "..", "examples", "oed.jl"))   # defines prob, fishing, measurements, oed, ps, st, sol, ...

using CairoMakie, LaTeXStrings
using SymbolicIndexingInterface: parameter_index, parameter_values_at_time

# States on a fine grid: the same layer with a finer saving grid and the optimized parameters.
_, reconstruct = pack(default_vectorizer(), ps)
fine = OEDLayer(remake(prob; saveat = 0.0:0.02:12.0), Tsit5(), fishing; params = [2, 3],
                measurements, controls = [1])
design, _ = fine(nothing, reconstruct(sol.u), st)
X = reduce(hcat, design.u)

# Piecewise constant controls on the grid of 48 intervals, read from the trajectory.
edges = collect(0.0:0.25:12.0)
value(s) = [parameter_values_at_time(design, design, t + 0.125)[parameter_index(design, s)]
            for t in edges[1:(end - 1)]]
step(v) = vcat(v, v[end])                     # stairs(:post) needs one value per edge

colors = Makie.wong_colors()
theme = merge(theme_latexfonts(), Theme(fontsize = 22))
with_theme(theme) do
    fig = Figure(; size = (1200, 420))
    states = Makie.Axis(fig[1, 1]; xlabel = L"t", ylabel = L"x(t)")
    for i in 1:2
        lines!(states, design.t, X[i, :]; color = colors[i], linewidth = 2)
    end
    controls = Makie.Axis(fig[1, 2]; xlabel = L"t", ylabel = L"u(t),\ w(t)")
    for (s, color, style) in ((:p₁, colors[3], :solid), (:w1, colors[1], :solid), (:w2, colors[2], :dash))
        stairs!(controls, edges, step(clamp.(value(s), 0, 1)); step = :post, color, linewidth = 2,
                linestyle = style)
    end
    Legend(fig[2, 1], [LineElement(; color = colors[i], linewidth = 2) for i in 1:2],
           [L"x_1", L"x_2"]; orientation = :horizontal, framevisible = false)
    Legend(fig[2, 2], [LineElement(; color = colors[3], linewidth = 2),
                       LineElement(; color = colors[1], linewidth = 2),
                       LineElement(; color = colors[2], linewidth = 2, linestyle = :dash)],
           [L"u", L"w_1", L"w_2"]; orientation = :horizontal, framevisible = false)
    out = joinpath(get(ENV, "FIGURE_DIR", mkpath(joinpath(@__DIR__, "..", "results"))),
                   "oed-optimal-design.pdf")
    save(out, fig)
    save(replace(out, ".pdf" => ".png"), fig)
    println("saved ", normpath(out))
end
