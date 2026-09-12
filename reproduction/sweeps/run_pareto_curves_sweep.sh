#!/usr/bin/env bash
# Pareto-curve scenarios (F1=Normal-vs-Pareto cosine, F2=Pareto-vs-Pareto cosine, F3=Normal-vs-Pareto piecewise)
# extended to cover all three methods (Wilcoxon, L_p, RKHS) at the
# framework's standard parameters (n=200, NREP=100, GAMMA=c(0,0.25,0.45),
# A_GRID, EC_GEOMETRY).
#
# Wilcoxon: 3 projections x m in {1,2,4} x 3 scenarios = 27 cells
# L_p     : 3 projections x p=2 x 3 scenarios          = 9 cells
# RKHS    : 3 projection-kernel pairs x 3 scenarios     = 9 cells
# Total: 45 driver invocations.
#
# Outputs land in `results/<method>_<scenario>_<tag>_<projection>[_m<m>].csv`.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/pareto_curves_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

run_driver () {
  local method=$1 scenario=$2 projection=$3 wkernel=${4:-h_w} kernel=${5:-gaussian} m=${6:-1} p=${7:-2}
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp method=$method scenario=$scenario projection=$projection wkernel=$wkernel kernel=$kernel m=$m p=$p ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- '$method'
    SCENARIO        <- '$scenario'
    PROJECTION      <- '$projection'
    WILCOXON_KERNEL <- '$wkernel'
    KERNEL          <- '$kernel'
    M               <- ${m}L
    P               <- ${p}L
    NCORES          <- 12L
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

# Wilcoxon x {ito, l2, c_max_var} x m in {1, 2, 4} x {F1, F2, F3}
for scenario in F1 F2 F3; do
  for proj in ito l2 c_max_var; do
    for m in 1 2 4; do
      run_driver wilcoxon "$scenario" "$proj" h_0 gaussian "$m"
    done
  done
done

# L_p x {ito, l2, c_max_var} x p=2 x {F1, F2, F3}
for scenario in F1 F2 F3; do
  for proj in ito l2 c_max_var; do
    run_driver lp "$scenario" "$proj"
  done
done

# RKHS x {ito/Gaussian, l2/min, c_max_var/Gaussian} x {F1, F2, F3}
for scenario in F1 F2 F3; do
  run_driver rkhs "$scenario" ito        h_w gaussian
  run_driver rkhs "$scenario" l2         h_w min
  run_driver rkhs "$scenario" c_max_var  h_w gaussian
done

echo "=== $(date '+%H:%M:%S') 45 Pareto-curve runs complete; CSVs in $ROOT/results/ ===" | tee -a "$LOG"
