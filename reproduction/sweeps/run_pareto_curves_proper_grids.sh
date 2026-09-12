#!/usr/bin/env bash
# Re-run Wilcoxon and L_p sweeps for F1/F2/F3 with per-scenario a-grids:
#   F1 (Normal-vs-Pareto, cosine):  a in {1.5, 2, 2.5, 3, 5, 10}    (a = 0 size check is run separately)
#   F2 (Pareto-vs-Pareto, cosine):  a in {0.5, 1, 1.5, 2, 5}        (a = 1.5 implicit size check)
#   F3 (Normal-vs-Pareto, piecewise): a in {0.5, 1, 1.5, 2, 5}
#
# Wilcoxon: 3 projections x min{1,2,4} x 3 scenarios = 27 cells
# L_p     : 3 projections x p=2     x 3 scenarios = 9 cells
#
# RKHS handled separately by run_pareto_curves_rkhs_calibrated.R (which
# also takes per-scenario a-grids).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/pareto_proper_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

A_F1='c(1.5, 2, 2.5, 3, 5, 10)'
A_F2='c(0.5, 1, 1.5, 2, 5)'
A_F3='c(0.5, 1, 1.5, 2, 2.5, 5)'

run () {
  local method=$1 scenario=$2 projection=$3 kernel=${4:-h_w} m=${5:-1} grid=$6
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp method=$method scenario=$scenario proj=$projection kernel=$kernel m=$m a=$grid ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- '$method'
    SCENARIO        <- '$scenario'
    PROJECTION      <- '$projection'
    WILCOXON_KERNEL <- '$kernel'
    KERNEL          <- 'gaussian'
    M               <- ${m}L
    A_GRID          <- $grid
    NCORES          <- 12L
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

# Wilcoxon x {ito, l2, c_max_var} x m in {1, 2, 4} x {F1, F2, F3}
for proj in ito l2 c_max_var; do
  for m in 1 2 4; do
    run wilcoxon F1 "$proj" h_0 "$m" "$A_F1"
    run wilcoxon F2 "$proj" h_0 "$m" "$A_F2"
    run wilcoxon F3 "$proj" h_0 "$m" "$A_F3"
  done
done

# L_p x {ito, l2, c_max_var} x p=2 x {F1, F2, F3}
for proj in ito l2 c_max_var; do
  run lp F1 "$proj" h_w 1 "$A_F1"
  run lp F2 "$proj" h_w 1 "$A_F2"
  run lp F3 "$proj" h_w 1 "$A_F3"
done

echo "=== $(date '+%H:%M:%S') Wilcoxon (27) + L_p (9) cells complete ===" | tee -a "$LOG"
echo "Now run: Rscript $ROOT/run_pareto_curves_rkhs_calibrated.R"
