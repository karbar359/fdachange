#!/usr/bin/env bash
# P4b2: extend the F4 L_p sweep to p in {4, 10} across the three projections.
# Six runs total. Mirrors the configuration / outputs of run_F4_sweep.sh
# (n=200, NREP=100, GAMMA_GRID=c(0,0.25,0.45), A_GRID, EC_GEOMETRY).
#
# Outputs land in `results/lp_F4_p<p>_<projection>.csv`.
#
# Usage:
#   bash reproduction/sweeps/run_F4_lp_extended.sh

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/F4_lp_extended_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

run_lp () {
  local proj=$1 p=$2
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp lp projection=$proj P=$p ===" | tee -a "$LOG"
  Rscript -e "
    METHOD     <- 'lp'
    SCENARIO   <- 'F4'
    PROJECTION <- '$proj'
    P          <- ${p}L
    NCORES     <- 12L
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

for proj in c_max_var ito l2; do
  for p in 4 10; do
    run_lp "$proj" "$p"
  done
done

echo "=== $(date '+%H:%M:%S') 6 lp pin{4,10} runs complete; CSVs in $ROOT/results/ ===" | tee -a "$LOG"
