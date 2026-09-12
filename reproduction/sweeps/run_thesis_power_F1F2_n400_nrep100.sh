#!/usr/bin/env bash
# Defence-deck F1/F2 rerun at n = 400, N_rep = 100: 18 cells total.
#
#   wilcoxon (driver.R):  F1,F2 x {ito, l2, c_max_var}, h_0, m = 1   (6 cells)
#   lp p=2   (driver.R):  F1,F2 x {ito, l2, c_max_var}              (6 cells)
#   rkhs (calibrated):    F1,F2 x {gaussian+ito, min+runif+l2,
#                                  gaussian+c_max_var}              (6 cells)
#
# Why: the printed thesis tables _pareto_F{1,2}_{wilcoxon,lp,rkhs}_*.tex are
# N_rep = 20; the thesis power figures use N_rep = 100 everywhere. The earlier
# wilcoxon-only rerun (run_pareto_F1F2_n400_nrep100.sh) covered m in {1,2,4}
# but its CSVs were deleted; this sweep recomputes the m = 1 cells (the m 2/4
# values are preserved in code/assemble_thesis_power.R) and adds lp/rkhs.
#
# RKHS goes through run_thesis_F1F2_rkhs_calibrated_nrep100.R, the
# rank-calibrated pipeline that produced the printed tables; the driver.R
# RKHS path is a different statistic and must not be used here.
#
# All outputs carry a thesis_ prefix so aggregate_pareto_curves_tex.R
# (glob ^(lp|wilcoxon|rkhs)_F[1-3]) can never pick them up when
# regenerating thesis tables. Idempotent: cells whose output CSV already
# exists and is non-empty are skipped.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/thesis_F1F2_nrep100_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG (REPRO_CACHE_DIR=$REPRO_CACHE_DIR)"

# Thesis (paper 1) a-grids; a is the Pareto shape in the change window.
A_F1='c(1.5, 2, 2.5, 3, 5, 10)'
A_F2='c(0.5, 1, 1.5, 2, 5)'
NCORES=16

run_driver () {
  local method=$1 sc=$2 proj=$3 a_grid=$4 tag=$5
  local out="results/thesis_${method}_${sc}_n400_nrep100_${tag}_${proj}.csv"
  if [ -s "$ROOT/$out" ]; then
    echo "SKIP $out (exists)" | tee -a "$LOG"
    return 0
  fi
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp $method $sc proj=$proj a=$a_grid -> $out ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- '$method'
    SCENARIO        <- '$sc'
    PROJECTION      <- '$proj'
    WILCOXON_KERNEL <- 'h_0'
    P               <- 2L
    M               <- 1L
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
    run_driver wilcoxon "$sc" "$proj" "$a_grid" h_0_beta0
  done
done

for sc in F1 F2; do
  a_grid=$A_F1; [ "$sc" = "F2" ] && a_grid=$A_F2
  for proj in ito l2 c_max_var; do
    run_driver lp "$sc" "$proj" "$a_grid" p2
  done
done

echo "=== $(date '+%H:%M:%S') driver cells done; starting calibrated RKHS ===" | tee -a "$LOG"
Rscript "$ROOT/run_thesis_F1F2_rkhs_calibrated_nrep100.R" 2>&1 | tee -a "$LOG"

echo "=== $(date '+%H:%M:%S') 18 cells complete ===" | tee -a "$LOG"
