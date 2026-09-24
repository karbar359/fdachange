#!/usr/bin/env bash
# Re-run the entire dissertation simulation suite at n=400 with NCORES=20.
# Resumable: cells whose output CSV already exists and is non-empty are
# skipped, so the script can be interrupted (Ctrl+C / SIGTERM) and re-run
# to pick up from approximately the same point.
#
# Coverage:
#  Pareto F1/F2/F3 :  Wilcoxon (3 proj x 3 m x 3 scen = 27)
#                  +  L_p      (3 proj x 3 scen      =  9)
#                  +  RKHS     (3 proj x 3 scen      =  9)  via R wrapper
#  F4              :  Wilcoxon (3 proj x 3 m         =  9)
#                  +  L_p      (3 proj               =  3)
#                  +  RKHS     (3 (proj,kernel)      =  3)
#                  +  Birbilas (1                    =  1)
#  Scalar S1..S5   :  Wilcoxon (5 scen x 3 m         = 15)
#                  +  L_p      (5 scen               =  5)
#                  +  RKHS     ((S1,gauss)+(S2,min)+(S4,2_max) = 3)
#  ----------------------------------------------------------------------
#  Total: 84 cells. Pre-flight: F4 RKHS finite null + Birbilas finite null
#  at n=400 (one-time, ~30 min).
#
# Outputs land in results/<method>_<scenario>_n400_<tag>.csv so the n=200
# baseline is untouched. The two finite-null caches at n=400 are appended
# to .raw/SD_h_F4_finite.csv and .raw/E_n_beta_gamma_finite.csv.
#
# Usage:
#   bash reproduction/sweeps/run_all_n400.sh
#   # background, survives terminal close:
#   nohup bash reproduction/sweeps/run_all_n400.sh >/dev/null 2>&1 &
#
# To resume after interruption: just re-run; everything done is skipped.

set -uo pipefail   # NOT -e: keep going through cell failures

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/all_n400_$(date '+%Y%m%d_%H%M%S').log"
NCORES=20
mkdir -p "$ROOT/results"

# Locate .raw/ (walk up from $ROOT until found)
RAW_DIR=""
cur="$ROOT"
for _ in 1 2 3 4 5 6; do
  if [ -d "$cur/.raw" ]; then RAW_DIR="$cur/.raw"; break; fi
  cur="$(dirname "$cur")"
done
if [ -z "$RAW_DIR" ]; then
  echo "ERROR: could not find .raw/ above $ROOT" >&2
  exit 1
fi

log () { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

log "Starting n=400 sweep. ROOT=$ROOT RAW=$RAW_DIR NCORES=$NCORES"
log "Log: $LOG"

# ---------------------------------------------------------------------------
# 0. Pre-flight: ensure finite-null caches at n=400 exist in .raw/
# ---------------------------------------------------------------------------

ensure_F4_RKHS_null () {
  local raw_csv="$RAW_DIR/SD_h_F4_finite.csv"
  if awk -F, 'NR>1 && $2=="400"{found=1; exit} END{exit !found}' "$raw_csv"; then
    log "F4 RKHS finite null at n=400 already in $raw_csv -- skip"
    return 0
  fi
  local out="$ROOT/results/SD_h_F4_finite_n400.csv"
  if [ ! -f "$out" ] || [ "$(wc -l < "$out")" -le 1 ]; then
    log "Generating F4 RKHS finite null at n=400..."
    local tmp_gen="$ROOT/results/.tmp_gen_F4_n400.R"
    sed -e 's/^N_S       <- c(200L)/N_S       <- c(400L)/' \
        -e "s|^OUT_CSV   <- NULL|OUT_CSV   <- '$out'|" \
        -e "s/^NCORES    <- max(.*)/NCORES    <- ${NCORES}L/" \
        "$ROOT/generators/rkhs_F4_null.R" > "$tmp_gen"
    # Place inside generators/ so script_path resolves root correctly.
    mv "$tmp_gen" "$ROOT/generators/.tmp_gen_F4_n400.R"
    Rscript "$ROOT/generators/.tmp_gen_F4_n400.R" 2>&1 | tee -a "$LOG"
    local rc=${PIPESTATUS[0]}
    rm -f "$ROOT/generators/.tmp_gen_F4_n400.R"
    if [ "$rc" -ne 0 ] || [ ! -f "$out" ]; then
      log "ERROR: F4 RKHS null generation failed (rc=$rc); aborting"
      exit 1
    fi
  else
    log "  found $out from a prior interrupted run"
  fi
  log "Appending n=400 rows to $raw_csv"
  awk 'NR>1' "$out" >> "$raw_csv"
  log "  done"
}

ensure_far1_null () {
  local raw_csv="$RAW_DIR/E_n_beta_gamma_finite.csv"
  if awk -F, 'NR>1 && $2=="400"{found=1; exit} END{exit !found}' "$raw_csv"; then
    log "Birbilas finite null at n=400 already in $raw_csv -- skip"
    return 0
  fi
  local out="$ROOT/results/E_n_beta_gamma_finite_n400.csv"
  if [ ! -f "$out" ] || [ "$(wc -l < "$out")" -le 1 ]; then
    log "Generating Birbilas finite null at n=400..."
    local tmp_gen="$ROOT/generators/.tmp_gen_far1_n400.R"
    sed -e 's/^N_S       <- c(100L, 200L, 500L)/N_S       <- c(400L)/' \
        -e "s|^OUT_CSV   <- NULL|OUT_CSV   <- '$out'|" \
        -e "s/^NCORES    <- max(.*)/NCORES    <- ${NCORES}L/" \
        "$ROOT/generators/far1_finite_null.R" > "$tmp_gen"
    Rscript "$tmp_gen" 2>&1 | tee -a "$LOG"
    local rc=${PIPESTATUS[0]}
    rm -f "$tmp_gen"
    if [ "$rc" -ne 0 ] || [ ! -f "$out" ]; then
      log "ERROR: Birbilas null generation failed (rc=$rc); aborting"
      exit 1
    fi
  else
    log "  found $out from a prior interrupted run"
  fi
  log "Appending n=400 rows to $raw_csv"
  awk 'NR>1' "$out" >> "$raw_csv"
  log "  done"
}

ensure_F4_RKHS_null
ensure_far1_null

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# True (zero) if cell needs to run; false (one) if output already exists.
needs_run () {
  local out="$1"
  if [ -f "$out" ] && [ "$(wc -l < "$out")" -gt 1 ]; then
    log "SKIP $(basename "$out") (exists)"
    return 1
  fi
  return 0
}

run_driver () {
  # Call driver.R via Rscript -e with pre-set globals. Uses .def() guards
  # so any pre-set variable wins over the driver default.
  # Args: out_csv, then any number of NAME=VALUE R-assignment strings.
  local out_csv="$1"; shift
  needs_run "$out_csv" || return 0
  local pre=""
  for kv in "$@"; do
    pre+="$kv; "
  done
  log "RUN $(basename "$out_csv")"
  Rscript -e "
    ${pre}
    NCORES <- ${NCORES}L
    OUT_CSV <- '$out_csv'
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

# ---------------------------------------------------------------------------
# 1. Pareto F1/F2/F3 -- Wilcoxon (27 cells)
# ---------------------------------------------------------------------------

A_F1='c(1.5, 2, 2.5, 3, 5, 10)'
A_F2='c(0.5, 1, 1.5, 2, 5)'
A_F3='c(0.5, 1, 1.5, 2, 2.5, 5)'
a_grid_for () {
  case "$1" in F1) echo "$A_F1" ;; F2) echo "$A_F2" ;; F3) echo "$A_F3" ;; esac
}

log "=== Pareto Wilcoxon (27 cells) ==="
for proj in ito l2 c_max_var; do
  for m in 1 2 4; do
    for scen in F1 F2 F3; do
      m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
      out="$ROOT/results/wilcoxon_${scen}_n400_h_0_beta0_${proj}${m_suffix}.csv"
      run_driver "$out" \
        "METHOD <- 'wilcoxon'" \
        "SCENARIO <- '$scen'" \
        "PROJECTION <- '$proj'" \
        "WILCOXON_KERNEL <- 'h_0'" \
        "M <- ${m}L" \
        "N <- 400L" \
        "A_GRID <- $(a_grid_for "$scen")"
    done
  done
done

# ---------------------------------------------------------------------------
# 2. Pareto F1/F2/F3 -- L_p (9 cells)
# ---------------------------------------------------------------------------

log "=== Pareto L_p (9 cells) ==="
for proj in ito l2 c_max_var; do
  for scen in F1 F2 F3; do
    out="$ROOT/results/lp_${scen}_n400_p2_${proj}.csv"
    run_driver "$out" \
      "METHOD <- 'lp'" \
      "SCENARIO <- '$scen'" \
      "PROJECTION <- '$proj'" \
      "P <- 2L" \
      "N <- 400L" \
      "A_GRID <- $(a_grid_for "$scen")"
  done
done

# ---------------------------------------------------------------------------
# 3. Pareto F1/F2/F3 -- RKHS (9 cells via dedicated R wrapper, idempotent)
# ---------------------------------------------------------------------------

log "=== Pareto RKHS calibrated (9 cells via R wrapper) ==="
Rscript "$ROOT/run_pareto_curves_rkhs_calibrated_n400.R" 2>&1 | tee -a "$LOG"

# ---------------------------------------------------------------------------
# 4. F4 cross-method (16 cells)
# ---------------------------------------------------------------------------

log "=== F4 (16 cells) ==="

# F4 L_p x 3 projections
for proj in c_max_var ito l2; do
  out="$ROOT/results/lp_F4_n400_p2_${proj}.csv"
  run_driver "$out" \
    "METHOD <- 'lp'" "SCENARIO <- 'F4'" \
    "PROJECTION <- '$proj'" "P <- 2L" "N <- 400L"
done

# F4 Wilcoxon x 3 projections x m in {1,2,4}
for proj in c_max_var ito l2; do
  for m in 1 2 4; do
    m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
    out="$ROOT/results/wilcoxon_F4_n400_h_0_beta0_${proj}${m_suffix}.csv"
    run_driver "$out" \
      "METHOD <- 'wilcoxon'" "SCENARIO <- 'F4'" \
      "PROJECTION <- '$proj'" "WILCOXON_KERNEL <- 'h_0'" \
      "M <- ${m}L" "N <- 400L"
  done
done

# F4 RKHS x 3 (projection, kernel) pairs
# F4 dispatches via load_rkhs_F4_null keyed on (n, gamma, projection, kernel)
# so the n=400 finite-null cache (just generated above) drives the cv.
run_driver "$ROOT/results/rkhs_F4_n400_gaussian_rnorm_c_max_var.csv" \
  "METHOD <- 'rkhs'" "SCENARIO <- 'F4'" \
  "PROJECTION <- 'c_max_var'" "KERNEL <- 'gaussian'" "N <- 400L"
run_driver "$ROOT/results/rkhs_F4_n400_gaussian_rnorm_ito.csv" \
  "METHOD <- 'rkhs'" "SCENARIO <- 'F4'" \
  "PROJECTION <- 'ito'" "KERNEL <- 'gaussian'" "N <- 400L"
run_driver "$ROOT/results/rkhs_F4_n400_min_rnorm_l2.csv" \
  "METHOD <- 'rkhs'" "SCENARIO <- 'F4'" \
  "PROJECTION <- 'l2'" "KERNEL <- 'min'" "N <- 400L"

# F4 Birbilas
run_driver "$ROOT/results/far1_F4_n400_beta1_cv_adjusted.csv" \
  "METHOD <- 'far1'" "SCENARIO <- 'F4'" \
  "PROJECTION <- 'none'" "FAR1_BETA_ID <- 'beta1'" \
  "FAR1_CV <- 'adjusted'" "N <- 400L"

# ---------------------------------------------------------------------------
# 5. Scalar S1-S5 (23 cells)
# ---------------------------------------------------------------------------

log "=== Scalar S1-S5: Wilcoxon (15 cells) ==="
for scen in S1 S2 S3 S4 S5; do
  for m in 1 2 4; do
    m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
    out="$ROOT/results/wilcoxon_${scen}_n400_h_w_beta0${m_suffix}.csv"
    run_driver "$out" \
      "METHOD <- 'wilcoxon'" "SCENARIO <- '$scen'" \
      "PROJECTION <- 'none'" "WILCOXON_KERNEL <- 'h_w'" \
      "M <- ${m}L" "N <- 400L"
  done
done

log "=== Scalar S1-S5: L_p (5 cells) ==="
for scen in S1 S2 S3 S4 S5; do
  out="$ROOT/results/lp_${scen}_n400_p2.csv"
  run_driver "$out" \
    "METHOD <- 'lp'" "SCENARIO <- '$scen'" \
    "PROJECTION <- 'none'" "P <- 2L" "N <- 400L"
done

log "=== Scalar S1/S2/S4: RKHS (3 cells) ==="
# Match sec_scalar.tex tab:pow-rkhs-all coverage: kernel paired with H_0 marginal.
#   S1 (Normal)         + gaussian   + baseline ""        (N(0,1) cache)
#   S2 (chi^2_1)        + min        + baseline rnorm.sq  (chi^2_1 cache)
#   S4 (Beta(1,9))      + 2_max      + baseline rbeta.1.9 (Beta(1,9) cache, calibrated)
run_driver "$ROOT/results/rkhs_S1_n400_gaussian_rnorm.csv" \
  "METHOD <- 'rkhs'" "SCENARIO <- 'S1'" \
  "PROJECTION <- 'none'" "KERNEL <- 'gaussian'" "BASELINE <- ''" "N <- 400L"
run_driver "$ROOT/results/rkhs_S2_n400_min_rnorm.sq.csv" \
  "METHOD <- 'rkhs'" "SCENARIO <- 'S2'" \
  "PROJECTION <- 'none'" "KERNEL <- 'min'" "BASELINE <- 'rnorm.sq'" "N <- 400L"
run_driver "$ROOT/results/rkhs_S4_n400_2_max_rbeta.1.9.csv" \
  "METHOD <- 'rkhs'" "SCENARIO <- 'S4'" \
  "PROJECTION <- 'none'" "KERNEL <- '2_max'" "BASELINE <- 'rbeta.1.9'" "N <- 400L"

log "=== ALL n=400 CELLS ATTEMPTED ==="
log "  Re-run this script to retry any cells that errored or were interrupted."
