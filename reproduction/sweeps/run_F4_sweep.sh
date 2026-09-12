#!/usr/bin/env bash
# Run all F4 cross-method paper-grade comparisons in one shot. Each
# invocation calls `driver.R` via `Rscript -e`, pre-setting METHOD /
# PROJECTION / WILCOXON_KERNEL in the global env so the .def guards in
# driver.R use them instead of the defaults.
#
# Configurations (16 runs total):
#   - lp x {c_max_var, ito, l2}        = 3 runs (P = 2, ecdf null, M=1)
#   - wilcoxon x {c_max_var, ito, l2}
#       x M in {1, 2, 4}                = 9 runs (kernel h_0; paper 1's
#     T^(d) outer-max aggregation across multiple chunk sizes.
#     N/M in {200, 100, 50} for n=200.)
#   - rkhs x {(c_max_var, gaussian),
#             (ito, gaussian),
#             (l2, min)}                = 3 runs (F4 finite-null cache,
#     M=1; gaussian for real-valued projections, min for R^+-valued l2)
#   - birbilas x none                   = 1 run (matrix-input, finite-null
#     cache, M=1)
#   - wiener / pc1 / pc_d projections   skipped (not in the published
#     cross-method comparison)
#
# Estimated runtime: ~10-20 min on 10 cores at n=200, NREP=100.
#
# Usage:
#   bash reproduction/sweeps/run_F4_sweep.sh
#
# Outputs land in `results/<method>_F4_<tag>_<projection>.csv` (one per
# run, never overwritten).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/F4_sweep_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

run_driver () {
  local method=$1 projection=$2 wkernel=${3:-h_w} kernel=${4:-gaussian} m=${5:-1}
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp method=$method projection=$projection wkernel=$wkernel kernel=$kernel m=$m ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- '$method'
    SCENARIO        <- 'F4'
    PROJECTION      <- '$projection'
    WILCOXON_KERNEL <- '$wkernel'
    KERNEL          <- '$kernel'
    M               <- ${m}L
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

# Lp x projections (kernel-free, M=1)
for proj in c_max_var ito l2; do
  run_driver lp "$proj"
done

# Wilcoxon x projections x M in {1, 2, 4} (paper 1 functional setup:
# kernel h_0; T^(d) outer-max aggregation at multiple chunk sizes).
for proj in c_max_var ito l2; do
  for m in 1 2 4; do
    run_driver wilcoxon "$proj" h_0 gaussian "$m"
  done
done

# RKHS x F4: hand-picked (projection, kernel) pairs, M=1.
run_driver rkhs c_max_var h_w gaussian
run_driver rkhs ito       h_w gaussian
run_driver rkhs l2        h_w min

# Birbilas: full-matrix statistic, no projection, M=1
run_driver birbilas none

echo "=== $(date '+%H:%M:%S') all 16 runs complete; CSVs in $ROOT/results/ ===" | tee -a "$LOG"
