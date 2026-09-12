#!/usr/bin/env bash
# Resume the Pareto-curve sweep for the missing RKHS cells. The first sweep
# halted at RKHS / L_2 / min kernel because the n=2000 cache doesn't have
# a "no-baseline" `min` column (only .runif and .rnorm.sq). Switch to the
# Gaussian kernel for all three projections, which the cache supports
# universally.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/pareto_curves_rkhs_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

run_rkhs () {
  local scenario=$1 projection=$2
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp rkhs scenario=$scenario projection=$projection kernel=gaussian ===" | tee -a "$LOG"
  Rscript -e "
    METHOD     <- 'rkhs'
    SCENARIO   <- '$scenario'
    PROJECTION <- '$projection'
    KERNEL     <- 'gaussian'
    M          <- 1L
    NCORES     <- 12L
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

# F1 ito completed in the prior run. Re-run F1 l2 + c_max_var, then F2x3, F3x3.
run_rkhs F1 l2
run_rkhs F1 c_max_var
for scenario in F2 F3; do
  for proj in ito l2 c_max_var; do
    run_rkhs "$scenario" "$proj"
  done
done

echo "=== $(date '+%H:%M:%S') 8 RKHS resume runs complete ===" | tee -a "$LOG"
