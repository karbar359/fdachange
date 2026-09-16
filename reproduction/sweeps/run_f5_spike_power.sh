#!/usr/bin/env bash
# F5 (Wiener + spike) power study: the experiment behind the thesis claim
# that point-evaluation projections suit spike-type (time-localised)
# changes. The spike is the dissertation's f_alpha function
# (1 - |t - 1/2|^alpha)^(1/alpha) at the printed alpha = 0.45 (the
# generate_f5 default). Sweeps the corrected Wilcoxon (h_0) and Lp
# (p = 2, ecdf) scans over the three thesis projections
# {c_max_var, l2, ito} at n = 400, N_rep = 100, driver-standard EC
# geometry (short l = 6 and long l = 80 centred windows) and a in
# {0, 0.5, 1, 2, 4, 10} (the archived smoke run used {0, 1, 10}). RKHS is
# excluded: no calibrated null for the projected F5 marginals.
#
# Aggregate the outputs with aggregators/aggregate_f5_spike.R.
# Existing outputs are skipped, so the script resumes after an interrupt.
#
# Usage:  bash reproduction/sweeps/run_f5_spike_power.sh

set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/f5_spike_power_$(date '+%Y%m%d_%H%M%S').log"
NCORES=12
mkdir -p "$ROOT/results"

log () { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

needs_run () {
  local f="$1"
  if [ -f "$f" ] && [ "$(wc -l < "$f")" -gt 1 ]; then
    log "SKIP $(basename "$f") (exists)"
    return 1
  fi
  return 0
}

run_driver () {
  local out_csv="$1"; shift
  needs_run "$out_csv" || return 0
  local pre=""
  for kv in "$@"; do pre+="$kv; "; done
  log "RUN $(basename "$out_csv")"
  Rscript -e "
    ${pre}
    EC_EXACT <- TRUE
    NCORES <- ${NCORES}L
    A_GRID <- c(0, 0.5, 1, 2, 4, 10)
    OUT_CSV <- '$out_csv'
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

log "Starting F5 spike power study. NCORES=$NCORES"

for proj in c_max_var l2 ito; do
  out="$ROOT/results/wilcoxon_F5_n400_h_0_beta0_corrected_${proj}.csv"
  run_driver "$out" \
    "METHOD <- 'wilcoxon'" "SCENARIO <- 'F5'" \
    "PROJECTION <- '$proj'" "WILCOXON_KERNEL <- 'h_0'" \
    "WILCOXON_SCAN <- 'corrected'" \
    "N <- 400L" "N_REP <- 100L"
done

for proj in c_max_var l2 ito; do
  out="$ROOT/results/lp_F5_n400_p2_corrected_${proj}.csv"
  run_driver "$out" \
    "METHOD <- 'lp'" "SCENARIO <- 'F5'" \
    "PROJECTION <- '$proj'" "P <- 2L" "LP_VARIANT <- 'corrected'" \
    "N <- 400L" "N_REP <- 100L"
done

log "Sweep done. Aggregate with: Rscript $ROOT/aggregators/aggregate_f5_spike.R"
log "Log: $LOG"
