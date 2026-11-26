# Corleone.jl: examples of the SoftwareX paper

Reproduces Section 3 and Appendix A of the paper *Corleone.jl: A modular shooting framework for
dynamic optimization and scientific machine learning in Julia*:

- optimal control of the Lotka–Volterra fishing problem with single shooting, multiple shooting,
  and integer controls (Section 3.1),
- a universal differential equation fitted with CorleoneInverse.jl (Section 3.2),
- optimal experimental design with CorleoneOED.jl (Section 3.3),

including the BenchmarkTools timings of Appendix A and the figures of Section 3. The results of
the paper are part of the repository in `results/`.
Corleone.jl itself is developed at <https://github.com/SciML/Corleone.jl>. This repository is archived
on Zenodo: [doi:10.5281/zenodo.24551106](https://doi.org/10.5281/zenodo.24551106) (all versions; the paper uses
v1.0.1, [doi:10.5281/zenodo.24553574](https://doi.org/10.5281/zenodo.24553574)).

## Layout

```
Project.toml, Manifest.toml   Julia environment (exact versions of the paper)
examples/                     one script per example of Section 3
figures/                      figures of Section 3
scripts/                      run.sh (everything), run_all.jl (all examples), scaling.sh
results/                      results of the paper (overwritten by the scripts)
```

## Running

Requires Julia 1.12.7 (e.g. `juliaup add 1.12.7`). Only registered packages of
the General registry are used; `Project.toml` and `Manifest.toml` pin the exact versions of the
paper (Appendix A), among them Corleone 2.1.0, CorleoneOED 2.1.0, and CorleoneInverse 0.1.0
(only CEnum differs: 0.5.0 instead of 1.0.0, required by the plotting packages).

```bash
scripts/run.sh
```

instantiates the environment, runs all examples, writes the timings, and creates the figures;
a run takes about 25 minutes after installation, most of it for the integer-control example.
Everything is written to `results/`, which holds the output of the run reported in the paper:

| File | Content |
| --- | --- |
| `benchmarks.log` | output of all examples (objective values, iterations, BenchmarkTools trials), `InteractiveUtils.versioninfo()`, and the summary table of Appendix A |
| `benchmarks.json` | BenchmarkTools trials (read with `BenchmarkTools.load`) |
| `scaling_threads_<n>.log` | scaling study of Appendix A with n threads (`scripts/scaling.sh`) |
| `fishing-single-shooting.pdf`, `fishing-multiple-shooting-initial.pdf` (and `.png`) | figures of Section 3.1 (log: `fishing_figure.log`) |
| `ude-fit.pdf`, `ude-fit.png` | figure of Section 3.2 (log: `ude_figure.log`) |
| `oed-optimal-design.pdf`, `oed-optimal-design.png` | figure of Section 3.3 (log: `oed_figure.log`) |

Objective values, iteration counts, and residuals are deterministic and agree with the paper;
computation times depend on the machine (the paper reports an AMD Ryzen 7 8700F with four
threads). The environment variables `RUN_MINLP=false` (skip the integer-control example),
`BENCHMARK_SAMPLES=n` (timing samples of that example, default 10), and `THREADS=n` (default 4)
shorten or change the run.

Single examples, from the root of the repository:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/fishing_single_shooting.jl
julia --project=. -t 4 examples/fishing_multiple_shooting.jl
julia --project=. examples/ude.jl
julia --project=. examples/oed.jl
julia --project=. -t 4 scripts/run_all.jl   # all examples (RUN_MINLP=true for the integer example)
scripts/scaling.sh    # scaling study of Appendix A with 1, 2, 4, 8 threads (THREADS="1 2")
julia --project=. -t 4 figures/fishing_figure.jl
julia --project=. figures/ude_figure.jl
julia --project=. figures/oed_figure.jl
```

With Docker (environment installed and precompiled at build time):

```bash
docker build -t corleone-paper .
docker run --rm -v "$PWD/results":/work/results corleone-paper
```

## Files

| File | Paper |
| --- | --- |
| `examples/fishing_single_shooting.jl` | Section 3.1, Listing 1 |
| `examples/fishing_multiple_shooting.jl` | Section 3.1, multiple shooting (four threads) |
| `examples/fishing_multiple_shooting_sparse_jacobian.jl` | Section 3.1, multiple shooting with sparse constraint Jacobian |
| `examples/fishing_multiple_shooting_sparse_hessian.jl` | Section 3.1, multiple shooting with sparse Jacobian and exact sparse Lagrangian Hessian |
| `examples/fishing_multiple_shooting_scaling.jl`, `scripts/scaling.sh` | scaling study of Appendix A (shooting nodes and threads; not part of `scripts/run.sh`) |
| `examples/fishing_minlp.jl` | Section 3.1, integer controls (Juniper and Ipopt) |
| `examples/fishing_cia.jl` | Section 3.1, sum-up rounding of the relaxed fishing control |
| `examples/ude.jl` | Section 3.2, Listing 2 |
| `examples/oed.jl` | Section 3.3, Listing 3 |
| `examples/oed_cia.jl` | Appendix A, sum-up rounding of the relaxed measurement weights |
| `scripts/run_all.jl` | runs the nine examples above and prints the summary table |
| `figures/fishing_figure.jl`, `figures/ude_figure.jl`, `figures/oed_figure.jl` | figures of Sections 3.1, 3.2, and 3.3 |
| `scripts/run.sh` | runs everything |

Each example script states the problem as in the paper's listing, solves it once, prints the key
numbers, and reports the wall-clock time of the solve with `BenchmarkTools.@benchmark`
(compilation and problem construction excluded). Differences to the listings in the paper:

- `ComponentArrays` is loaded in every script, because Corleone uses its extension for the
  default parameter vectorizer of `OptimizationProblem`;
- Ipopt runs with `tol = 1e-6` and `print_level = 0`, and Juniper with its log switched off;
- `oed.jl` also computes the uniform baseline design and the improvement factor;
- `ude.jl` includes the data generation and the scaling of the initial weights.

## Citation

See `CITATION.cff`.

## License

MIT, see `LICENSE`.
