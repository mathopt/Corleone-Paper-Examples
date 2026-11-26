#!/usr/bin/env bash
# Runs all examples of Section 3 of the Corleone.jl paper with BenchmarkTools timings and creates
# the figures of Section 3. All output goes to RESULTS_DIR (default: results/).
#   RUN_MINLP=false        skip the integer-control example (about 1.5 min per solve)
#   BENCHMARK_SAMPLES=n    timing samples of the integer-control example (default 10)
#   THREADS=n              Julia threads of the examples (default 4, as in the paper)
# Run:  scripts/run.sh
set -euo pipefail
cd "$(dirname "$0")/.."
export RESULTS_DIR="${RESULTS_DIR:-$PWD/results}"
export FIGURE_DIR="${FIGURE_DIR:-$RESULTS_DIR}"
mkdir -p "$RESULTS_DIR"

julia --version
julia --project=. -e 'using Pkg; Pkg.instantiate()'

RUN_MINLP="${RUN_MINLP:-true}" BENCHMARK_SAMPLES="${BENCHMARK_SAMPLES:-10}" \
    julia --project=. --threads="${THREADS:-4}" scripts/run_all.jl 2>&1 | tee "$RESULTS_DIR/benchmarks.log"
julia --project=. --threads="${THREADS:-4}" figures/fishing_figure.jl 2>&1 | tee "$RESULTS_DIR/fishing_figure.log"
julia --project=. figures/ude_figure.jl 2>&1 | tee "$RESULTS_DIR/ude_figure.log"
julia --project=. figures/oed_figure.jl 2>&1 | tee "$RESULTS_DIR/oed_figure.log"
