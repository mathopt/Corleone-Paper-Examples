#!/usr/bin/env bash
# Runs examples/fishing_multiple_shooting_scaling.jl with 1, 2, 4, and 8 threads (override with
# THREADS) and writes one log per thread count to RESULTS_DIR (default: results/).
# Run:  scripts/scaling.sh
set -euo pipefail
cd "$(dirname "$0")/.."
out="${RESULTS_DIR:-results}"
mkdir -p "$out"
for t in ${THREADS:-1 2 4 8}; do
    julia --project="${JULIA_PROJECT:-.}" --threads="$t" examples/fishing_multiple_shooting_scaling.jl \
        2>&1 | tee "$out/scaling_threads_$t.log"
done
