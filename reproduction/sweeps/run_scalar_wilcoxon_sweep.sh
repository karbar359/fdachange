#!/usr/bin/env bash
# P4b3: scalar Wilcoxon sweep over S1..S5 x min{1,2,4}.
# Mirrors the L_p / RKHS scalar evaluation grid (n=200, NREP=100,
# GAMMA_GRID = c(0,0.25,0.45), A_GRID, EC_GEOMETRY) so the new table
# slots into sec_scalar.tex alongside the existing L_p and RKHS tables.
#
# Outputs land in `results/wilcoxon_<scenario>_h_w_beta0[_m<m>].csv`.
# Total 15 runs.
#
# Usage:
#   bash reproduction/sweeps/run_scalar_wilcoxon_sweep.sh

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/scalar_wilcoxon_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

run_w () {
  local scenario=$1 m=$2
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp wilcoxon scenario=$scenario m=$m ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- 'wilcoxon'
    SCENARIO        <- '$scenario'
    PROJECTION      <- 'none'
    WILCOXON_KERNEL <- 'h_w'
    M               <- ${m}L
    NCORES          <- 12L
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

for scenario in S1 S2 S3 S4 S5; do
  for m in 1 2 4; do
    run_w "$scenario" "$m"
  done
done

echo "=== $(date '+%H:%M:%S') 15 scalar Wilcoxon runs complete; CSVs in $ROOT/results/ ===" | tee -a "$LOG"
