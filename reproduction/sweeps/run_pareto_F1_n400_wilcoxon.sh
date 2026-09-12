#!/usr/bin/env bash
# F1 (Normal-vs-Pareto, cosine) Wilcoxon at n=400, all 3 projections x min{1,2,4} = 9 cells.
# Outputs go to results/wilcoxon_F1_n400_<tag>.csv so they don't clobber
# the n=200 baseline.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/pareto_F1_n400_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

run_w () {
  local proj=$1 m=$2
  local suffix=""; [ "$m" -gt 1 ] && suffix="_m${m}"
  local out="results/wilcoxon_F1_n400_h_0_beta0_${proj}${suffix}.csv"
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp wilcoxon F1 proj=$proj m=$m -> $out ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- 'wilcoxon'
    SCENARIO        <- 'F1'
    PROJECTION      <- '$proj'
    WILCOXON_KERNEL <- 'h_0'
    M               <- ${m}L
    N               <- 400L
    NCORES          <- 12L
    OUT_CSV         <- '$ROOT/$out'
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

for proj in ito l2 c_max_var; do
  for m in 1 2 4; do
    run_w "$proj" "$m"
  done
done

echo "=== $(date '+%H:%M:%S') 9 cells complete ===" | tee -a "$LOG"
