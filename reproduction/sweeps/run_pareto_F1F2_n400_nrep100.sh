#!/usr/bin/env bash
# F1 and F2 (cosine Pareto scenarios) Wilcoxon at n=400 with N_rep=100,
# all 3 projections x m in {1,2,4} = 18 cells. The thesis tables for F1/F2
# were built from N_rep=20 runs; this sweep recomputes them at N_rep=100
# for the thesis power plots. Outputs carry an explicit _nrep100 tag so the
# N_rep=20 baseline CSVs are not clobbered.
#
# A_GRID must be set per scenario: the driver's default is the scalar
# (Ch.4/Ch.5) grid {0, 0.25, 0.5, 1, 2, 4, 8}, but the thesis F1/F2 tables
# sweep the paper-1 grids below. Without the override the sweep runs on the
# wrong a-values and cannot be compared against the thesis tables.
#
# Null caches are read from .raw/ via REPRO_CACHE_DIR; when run from a worktree
# that has no .raw/, point it at the main checkout (default below).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/pareto_F1F2_n400_nrep100_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG (REPRO_CACHE_DIR=$REPRO_CACHE_DIR)"

# Thesis (paper 1) a-grids; a is the Pareto shape in the change window.
A_F1='c(1.5, 2, 2.5, 3, 5, 10)'
A_F2='c(0.5, 1, 1.5, 2, 5)'
NCORES=16

run_w () {
  local sc=$1 proj=$2 m=$3 a_grid=$4
  local suffix=""; [ "$m" -gt 1 ] && suffix="_m${m}"
  local out="results/wilcoxon_${sc}_n400_nrep100_h_0_beta0_${proj}${suffix}.csv"
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp wilcoxon $sc proj=$proj m=$m a=$a_grid -> $out ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- 'wilcoxon'
    SCENARIO        <- '$sc'
    PROJECTION      <- '$proj'
    WILCOXON_KERNEL <- 'h_0'
    M               <- ${m}L
    A_GRID          <- $a_grid
    N               <- 400L
    N_REP           <- 100L
    NCORES          <- ${NCORES}L
    OUT_CSV         <- '$ROOT/$out'
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

for sc in F1 F2; do
  a_grid=$A_F1; [ "$sc" = "F2" ] && a_grid=$A_F2
  for proj in ito l2 c_max_var; do
    for m in 1 2 4; do
      run_w "$sc" "$proj" "$m" "$a_grid"
    done
  done
done

echo "=== $(date '+%H:%M:%S') 18 cells complete ===" | tee -a "$LOG"
