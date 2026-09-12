#!/usr/bin/env bash
# Re-run every Wilcoxon cell of the simulation chapter with the corrected scan
# (WILCOXON_SCAN='exact'), matching the printed parameterisation cell for cell.
#
# The printed tables were produced by the legacy scan, which drops column 1
# from the complement and so understates the statistic by roughly 1 per cent at
# n = 400 (see methods/wilcoxon/run.R). Outputs here carry an `_exact` tag so
# they never share a filename with the legacy CSVs, and so that
# aggregate_pareto_curves_tex.R's glob cannot pick them up by accident.
#
# Parameterisation is taken from the CSVs that fed the print:
#   scalar S1-S5  h_w, m in {1,2,4}, N_rep=100, A = default
#   F1            h_0, 3 projections x m in {1,2,4}, N_rep=20,  A = 1.5..10
#   F2            h_0, 3 projections x m in {1,2,4}, N_rep=20,  A = 0.5..5
#   F3            h_0, 3 projections x m in {1,2,4}, N_rep=100, A = 0.5..5
#   F4            h_0, 3 projections x m in {1,2,4}, N_rep=100, A = default
#
# Existing outputs are skipped, so the script can be re-run after an interrupt.
#
# Usage:  bash reproduction/sweeps/run_exact_scan_n400.sh

set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/exact_scan_n400_$(date '+%Y%m%d_%H%M%S').log"
NCORES=20
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
    WILCOXON_SCAN <- 'exact'
    NCORES <- ${NCORES}L
    OUT_CSV <- '$out_csv'
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

log "Starting exact-scan Wilcoxon sweep. NCORES=$NCORES"

A_F1='c(1.5, 2, 2.5, 3, 5, 10)'
A_F2='c(0.5, 1, 1.5, 2, 5)'
A_F3='c(0.5, 1, 1.5, 2, 2.5, 5)'

# ---------------------------------------------------------------------------
# 1. Functional F1/F2/F3 (27 cells)
# ---------------------------------------------------------------------------
log "=== Functional F1/F2/F3 Wilcoxon, exact scan (27 cells) ==="
for proj in ito l2 c_max_var; do
  for m in 1 2 4; do
    for scen in F1 F2 F3; do
      m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
      case "$scen" in
        F1) a_grid="$A_F1"; nrep=20  ;;
        F2) a_grid="$A_F2"; nrep=20  ;;
        F3) a_grid="$A_F3"; nrep=100 ;;
      esac
      out="$ROOT/results/wilcoxon_${scen}_n400_h_0_beta0_exact_${proj}${m_suffix}.csv"
      run_driver "$out" \
        "METHOD <- 'wilcoxon'" \
        "SCENARIO <- '$scen'" \
        "PROJECTION <- '$proj'" \
        "WILCOXON_KERNEL <- 'h_0'" \
        "M <- ${m}L" \
        "N <- 400L" \
        "N_REP <- ${nrep}L" \
        "A_GRID <- $a_grid"
    done
  done
done

# ---------------------------------------------------------------------------
# 2. F4 (9 cells)
# ---------------------------------------------------------------------------
log "=== F4 Wilcoxon, exact scan (9 cells) ==="
for proj in ito l2 c_max_var; do
  for m in 1 2 4; do
    m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
    out="$ROOT/results/wilcoxon_F4_n400_h_0_beta0_exact_${proj}${m_suffix}.csv"
    run_driver "$out" \
      "METHOD <- 'wilcoxon'" "SCENARIO <- 'F4'" \
      "PROJECTION <- '$proj'" "WILCOXON_KERNEL <- 'h_0'" \
      "M <- ${m}L" "N <- 400L"
  done
done

# ---------------------------------------------------------------------------
# 3. Scalar S1-S5 (15 cells)
# ---------------------------------------------------------------------------
log "=== Scalar S1-S5 Wilcoxon, exact scan (15 cells) ==="
for scen in S1 S2 S3 S4 S5; do
  for m in 1 2 4; do
    m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
    out="$ROOT/results/wilcoxon_${scen}_n400_h_w_beta0_exact${m_suffix}.csv"
    run_driver "$out" \
      "METHOD <- 'wilcoxon'" "SCENARIO <- '$scen'" \
      "PROJECTION <- 'none'" "WILCOXON_KERNEL <- 'h_w'" \
      "M <- ${m}L" "N <- 400L"
  done
done

log "Done. Log: $LOG"
