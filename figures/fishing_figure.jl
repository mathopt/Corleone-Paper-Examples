# Figures of Section 3.1 (figures/fishing-single-shooting.pdf and
# figures/fishing-multiple-shooting-initial.pdf), half page width each. From top to bottom:
# states x₁, x₂, fishing control u, and accumulated cost x₃.
# (a) optimal solution of the single-shooting problem (examples/fishing_single_shooting.jl);
# (b) initialization of the multiple-shooting problem (examples/fishing_multiple_shooting.jl):
#     shaded areas mark the shooting blocks, dotted lines the continuity residuals.
# Each example runs in its own module, which also solves it and times the solve.
#
# Run:  julia +1.12 --project=. -t 4 figures/fishing_figure.jl   (output: FIGURE_DIR, default: results/)

using Corleone, OrdinaryDiffEqTsit5
using CairoMakie, LaTeXStrings
using SymbolicIndexingInterface: parameter_index, parameter_values_at_time

function load_example(name)
    m = Module(Symbol(:Figure_, name))
    Base.include(m, joinpath(@__DIR__, "..", "examples", name * ".jl"))
    return m
end
ss = load_example("fishing_single_shooting")
ms = load_example("fishing_multiple_shooting")

# Trajectory on a fine saving grid: the example's layer, called with its problem remade with
# `saveat` as input (same segments, parameters, and state).
fine(m, p, st) = first(m.layer(Corleone.SciMLBase.remake(m.prob; saveat = 0.0:0.02:12.0), p, st))

# (a) Optimal single-shooting solution.
_, reconstruct_ss = pack(default_vectorizer(), ss.dynprob.ps)
ss_traj = fine(ss, reconstruct_ss(ss.sol.u), ss.opt.p)

# (b) Initialization of multiple shooting: the lifted states keep their initial values.
ms_traj = fine(ms, ms.ps, ms.st)

# Split a trajectory into segments at the duplicated times of the shooting nodes.
function segments(traj)
    t = traj.t
    cuts = [0; findall(i -> t[i] == t[i + 1], 1:(length(t) - 1)); length(t)]
    return [(cuts[k] + 1):cuts[k + 1] for k in 1:(length(cuts) - 1)]
end

edges = collect(0.0:0.25:12.0)
control_values(traj) = [parameter_values_at_time(traj, traj, t + 0.125)[parameter_index(traj, :p₂)]
                        for t in edges[1:(end - 1)]]

colors = Makie.wong_colors()
function fishing_figure(traj; nodes = Float64[])
    X = reduce(hcat, traj.u)
    fig = Figure(; size = (600, 640))
    axes = [Makie.Axis(fig[i, 1]; xlabel = i == 3 ? L"t" : "", ylabel = label)
            for (i, label) in enumerate((L"x_1(t),\ x_2(t)", L"u(t)", L"x_3(t)"))]
    linkxaxes!(axes...)
    blocks = [0.0; nodes; 12.0]
    for ax in axes, k in 1:2:(length(blocks) - 1)
        length(nodes) > 0 && vspan!(ax, blocks[k], blocks[k + 1]; color = (:gray, 0.12))
    end
    for seg in segments(traj)
        for (ax, i) in ((axes[1], 1), (axes[1], 2), (axes[3], 3))
            lines!(ax, traj.t[seg], X[i, seg]; color = colors[i], linewidth = 2)
        end
    end
    for seg in segments(traj)[1:(end - 1)]       # continuity residuals at the shooting nodes
        j = last(seg)
        for (ax, i) in ((axes[1], 1), (axes[1], 2), (axes[3], 3))
            lines!(ax, [traj.t[j], traj.t[j]], [X[i, j], X[i, j + 1]]; color = colors[i],
                   linestyle = :dot, linewidth = 2)
        end
    end
    u = control_values(traj)
    stairs!(axes[2], edges, vcat(u, u[end]); step = :post, color = colors[4], linewidth = 2)
    ylims!(axes[2], -0.05, 1.05)
    Legend(fig[4, 1], [LineElement(; color = colors[i], linewidth = 2) for i in (1, 2, 4, 3)],
           [L"x_1", L"x_2", L"u", L"x_3"]; orientation = :horizontal, framevisible = false)
    return fig
end

theme = merge(theme_latexfonts(), Theme(fontsize = 22))
outdir = get(ENV, "FIGURE_DIR", mkpath(joinpath(@__DIR__, "..", "results")))
with_theme(theme) do
    for (name, fig) in (("fishing-single-shooting", fishing_figure(ss_traj)),
                        ("fishing-multiple-shooting-initial",
                         fishing_figure(ms_traj; nodes = ms.shooting_points)))
        out = joinpath(outdir, name * ".pdf")
        save(out, fig)
        save(replace(out, ".pdf" => ".png"), fig)
        println("saved ", normpath(out))
    end
end
