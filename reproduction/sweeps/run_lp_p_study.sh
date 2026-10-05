#!/usr/bin/env bash
# Lp exponent study: how does the choice of p in the Lp scan affect power
# on the scalar scenarios? The thesis printed p = 2 only (legacy family);
# this sweep runs the CORRECTED Lp scan at p in {2, 4, 10} on S1-S5 so the
# comparison is like for like (the shipped null caches cover all three p).
# n = 400, N_rep = 100, gammas {0, 0.25, 0.45}, ch.4/5 a-grid, both
# centred EC geometries (short l = floor(log n)+1, long l = n/5), seed 42.
#
# Aggregate the outputs with aggregators/aggregate_lp_p.R.
# Existing outputs are skipped, so the script resumes after an interrupt.
#
# Usage:  bash reproduction/sweeps/run_lp_p_study.sh

set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/lp_p_study_$(date '+%Y%m%d_%H%M%S').log"
NCORES=12
# Keep the machine usable while the sweep runs: workers at lowest CPU
# priority, and any threaded BLAS capped to one thread per worker.
NICE="nice -n 19"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
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
  $NICE Rscript -e "
    ${pre}
    EC_EXACT <- TRUE
    NCORES <- ${NCORES}L
    OUT_CSV <- '$out_csv'
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

log "Starting Lp exponent study. NCORES=$NCORES"

for sc in S1 S2 S3 S4 S5; do
  for p in 2 4 10; do
    out="$ROOT/results/lp_${sc}_n400_p${p}_corrected.csv"
    run_driver "$out" \
      "METHOD <- 'lp'" "SCENARIO <- '$sc'" \
      "P <- ${p}L" "LP_VARIANT <- 'corrected'" \
      "N <- 400L" "N_REP <- 100L"
  done
done

log "Sweep done. Aggregate with: Rscript $ROOT/aggregators/aggregate_lp_p.R"
log "Log: $LOG"
