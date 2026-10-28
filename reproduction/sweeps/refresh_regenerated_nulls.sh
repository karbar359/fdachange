#!/usr/bin/env bash
# Re-run exactly the driver cells whose null cache was regenerated on
# 2026-10-28 with independent draws (see SEEDS.md), so the committed tables
# under data/ can be patched row by row instead of re-running the whole
# simulation chapter:
#
#   PATCH=1 Rscript compare_exact_vs_legacy.R                (FAMILY exact)
#   PATCH=1 Rscript -e "FAMILY <- 'corrected'; source('compare_exact_vs_legacy.R')"
#   PATCH=1 Rscript aggregators/assemble_thesis_power.R
#
# Cells, mirroring run_all_n400.sh, run_exact_scan_n400.sh,
# run_corrected_n400.sh and run_thesis_power_F1F2_n400_nrep100.sh tag for
# tag (N_rep 20 for F1/F2, 100 elsewhere; EC_EXACT for the corrected family):
#   * Wilcoxon, legacy + exact + corrected families, GAMMA_GRID = 0 only
#     (the regenerated T2 gamma = 0 cell): F1-F3 x 3 projections x m 1,2,4,
#     F4 x 3 x m, S1-S5 x m                          -> 153 driver runs
#   * RKHS corrected, every scenario (SD_h_n2000_corrected.csv): S1, S2, S4
#     through the driver, F1-F3 through run_pareto_curves_rkhs_corrected_n400.R
#   * RKHS F4 legacy + corrected, far1 F4 legacy + corrected (the finite F4
#     caches SD_h_F4_finite*.csv and E_n_beta_gamma_finite*.csv)
#   * thesis F1/F2 Wilcoxon reruns (N_rep = 100, every gamma)
#   * stability study, three families, GAMMA_GRID = 0
#
# Restricting a run to gamma = 0 changes which L'Ecuyer substream each
# (gamma, ec, a) task receives, so those rows carry fresh MC noise; both
# sides of every pair are re-run here, so the paired deltas stay exact.
#
# Outputs: results/refresh_2026-10-28/ (driver cells, stability summaries)
# and results/thesis_*.csv (the assembler's fixed location). Existing
# non-empty outputs are skipped, so the script resumes after an interrupt.
#
# Usage:  NCORES=16 bash reproduction/sweeps/refresh_regenerated_nulls.sh
#         (about 2 hours on 16 cores; every Rscript runs under nice)

set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT/results/refresh_2026-10-28"
NCORES="${NCORES:-16}"
NICE="${NICE:-10}"
mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/refresh_$(date '+%Y%m%d_%H%M%S').log"

log () { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

needs_run () {
  local f="$1"
  if [ -f "$f" ] && [ "$(wc -l < "$f")" -gt 1 ]; then
    log "SKIP $(basename "$f") (exists)"
    return 1
  fi
  return 0
}

# run_driver OUT_CSV "R assignment" ...  (pre-set globals win over the
# driver's .def defaults).
run_driver () {
  local out_csv="$1"; shift
  needs_run "$out_csv" || return 0
  local pre=""
  for kv in "$@"; do pre+="$kv; "; done
  log "RUN $(basename "$out_csv")"
  nice -n "$NICE" Rscript -e "
    ${pre}
    NCORES <- ${NCORES}L
    OUT_CSV <- '$out_csv'
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

A_F1='c(1.5, 2, 2.5, 3, 5, 10)'
A_F2='c(0.5, 1, 1.5, 2, 5)'
A_F3='c(0.5, 1, 1.5, 2, 2.5, 5)'
a_grid_for () { case "$1" in F1) echo "$A_F1" ;; F2) echo "$A_F2" ;; F3) echo "$A_F3" ;; esac; }
nrep_for    () { case "$1" in F1|F2) echo 20 ;; *) echo 100 ;; esac; }

family_pre () {
  case "$1" in
    legacy)    echo "WILCOXON_SCAN <- 'legacy'" ;;
    exact)     echo "WILCOXON_SCAN <- 'exact'" ;;
    corrected) echo "WILCOXON_SCAN <- 'corrected'; EC_EXACT <- TRUE" ;;
  esac
}
family_tag () { case "$1" in legacy) echo "" ;; exact) echo "_exact" ;; corrected) echo "_corrected" ;; esac; }

log "Refreshing cells that read regenerated null caches. NCORES=$NCORES OUT_DIR=$OUT_DIR"

# ---------------------------------------------------------------------------
# 1. Wilcoxon, three families, gamma = 0 (51 runs per family)
# ---------------------------------------------------------------------------
for fam in legacy exact corrected; do
  pre="$(family_pre "$fam")"; ftag="$(family_tag "$fam")"
  log "=== Wilcoxon $fam, gamma = 0: F1-F3 (27), F4 (9), S1-S5 (15) ==="
  for proj in ito l2 c_max_var; do
    for m in 1 2 4; do
      m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
      for scen in F1 F2 F3; do
        run_driver "$OUT_DIR/wilcoxon_${scen}_n400_h_0_beta0${ftag}_${proj}${m_suffix}.csv" \
          "$pre" "GAMMA_GRID <- 0" \
          "METHOD <- 'wilcoxon'" "SCENARIO <- '$scen'" \
          "PROJECTION <- '$proj'" "WILCOXON_KERNEL <- 'h_0'" \
          "M <- ${m}L" "N <- 400L" "N_REP <- $(nrep_for "$scen")L" \
          "A_GRID <- $(a_grid_for "$scen")"
      done
      run_driver "$OUT_DIR/wilcoxon_F4_n400_h_0_beta0${ftag}_${proj}${m_suffix}.csv" \
        "$pre" "GAMMA_GRID <- 0" \
        "METHOD <- 'wilcoxon'" "SCENARIO <- 'F4'" \
        "PROJECTION <- '$proj'" "WILCOXON_KERNEL <- 'h_0'" \
        "M <- ${m}L" "N <- 400L"
    done
  done
  for scen in S1 S2 S3 S4 S5; do
    for m in 1 2 4; do
      m_suffix=""; [ "$m" -gt 1 ] && m_suffix="_m${m}"
      run_driver "$OUT_DIR/wilcoxon_${scen}_n400_h_w_beta0${ftag}${m_suffix}.csv" \
        "$pre" "GAMMA_GRID <- 0" \
        "METHOD <- 'wilcoxon'" "SCENARIO <- '$scen'" \
        "PROJECTION <- 'none'" "WILCOXON_KERNEL <- 'h_w'" \
        "M <- ${m}L" "N <- 400L"
    done
  done
done

# ---------------------------------------------------------------------------
# 2. RKHS corrected: S1/S2/S4 via the driver, F1-F3 via the Pareto script
# ---------------------------------------------------------------------------
log "=== RKHS corrected S1/S2/S4 (3) ==="
run_driver "$OUT_DIR/rkhs_S1_n400_gaussian_rnorm_corrected.csv" \
  "EC_EXACT <- TRUE" "METHOD <- 'rkhs'" "SCENARIO <- 'S1'" "PROJECTION <- 'none'" \
  "KERNEL <- 'gaussian'" "BASELINE <- ''" "RKHS_VARIANT <- 'corrected'" "N <- 400L"
run_driver "$OUT_DIR/rkhs_S2_n400_min_rnorm.sq_corrected.csv" \
  "EC_EXACT <- TRUE" "METHOD <- 'rkhs'" "SCENARIO <- 'S2'" "PROJECTION <- 'none'" \
  "KERNEL <- 'min'" "BASELINE <- 'rnorm.sq'" "RKHS_VARIANT <- 'corrected'" "N <- 400L"
run_driver "$OUT_DIR/rkhs_S4_n400_2_max_rbeta.1.9_corrected.csv" \
  "EC_EXACT <- TRUE" "METHOD <- 'rkhs'" "SCENARIO <- 'S4'" "PROJECTION <- 'none'" \
  "KERNEL <- '2_max'" "BASELINE <- 'rbeta.1.9'" "RKHS_VARIANT <- 'corrected'" "N <- 400L"

log "=== RKHS corrected F1-F3 (9) via run_pareto_curves_rkhs_corrected_n400.R ==="
nice -n "$NICE" Rscript -e "
  NCORES <- ${NCORES}L
  RES_DIR <- '$OUT_DIR'
  source('$ROOT/run_pareto_curves_rkhs_corrected_n400.R')
" 2>&1 | tee -a "$LOG"

# ---------------------------------------------------------------------------
# 3. F4: RKHS legacy + corrected (3 + 3), far1 legacy + corrected (1 + 1)
# ---------------------------------------------------------------------------
log "=== F4 RKHS and far1, both families (8) ==="
for proj in c_max_var ito l2; do
  kern="gaussian"; [ "$proj" = "l2" ] && kern="min"
  run_driver "$OUT_DIR/rkhs_F4_n400_${kern}_rnorm_${proj}.csv" \
    "METHOD <- 'rkhs'" "SCENARIO <- 'F4'" \
    "PROJECTION <- '$proj'" "KERNEL <- '$kern'" "N <- 400L"
  run_driver "$OUT_DIR/rkhs_F4_n400_${kern}_rnorm_corrected_${proj}.csv" \
    "EC_EXACT <- TRUE" "METHOD <- 'rkhs'" "SCENARIO <- 'F4'" \
    "PROJECTION <- '$proj'" "KERNEL <- '$kern'" \
    "RKHS_VARIANT <- 'corrected'" "N <- 400L"
done
run_driver "$OUT_DIR/far1_F4_n400_beta1_cv_adjusted.csv" \
  "METHOD <- 'far1'" "SCENARIO <- 'F4'" \
  "PROJECTION <- 'none'" "FAR1_BETA_ID <- 'beta1'" \
  "FAR1_CV <- 'adjusted'" "N <- 400L"
run_driver "$OUT_DIR/far1_F4_n400_beta1_cv_adjusted_corrected.csv" \
  "EC_EXACT <- TRUE" "METHOD <- 'far1'" "SCENARIO <- 'F4'" \
  "PROJECTION <- 'none'" "FAR1_BETA_ID <- 'beta1'" \
  "FAR1_CV <- 'adjusted'" "FAR1_VARIANT <- 'corrected'" "N <- 400L"

# ---------------------------------------------------------------------------
# 4. Thesis F1/F2 Wilcoxon reruns (6), into the assembler's fixed location
# ---------------------------------------------------------------------------
log "=== Thesis F1/F2 Wilcoxon reruns (6) ==="
for sc in F1 F2; do
  for proj in ito l2 c_max_var; do
    run_driver "$ROOT/results/thesis_wilcoxon_${sc}_n400_nrep100_h_0_beta0_${proj}.csv" \
      "METHOD <- 'wilcoxon'" "SCENARIO <- '$sc'" "PROJECTION <- '$proj'" \
      "WILCOXON_KERNEL <- 'h_0'" "P <- 2L" "M <- 1L" \
      "A_GRID <- $(a_grid_for "$sc")" "N <- 400L" "N_REP <- 100L"
  done
done

# ---------------------------------------------------------------------------
# 5. Stability study, three families, gamma = 0
# ---------------------------------------------------------------------------
for fam in legacy exact corrected; do
  suffix=""; [ "$fam" != "legacy" ] && suffix="_${fam}"
  out="$OUT_DIR/stability_summary${suffix}.csv"
  needs_run "$out" || continue
  log "=== Stability, $fam scan, gamma = 0 ==="
  ec_exact=""; [ "$fam" = "corrected" ] && ec_exact="EC_EXACT <- TRUE;"
  nice -n "$NICE" Rscript -e "
    WILCOXON_SCAN <- '$fam'; $ec_exact
    GAMMA_GRID <- 0; NCORES <- ${NCORES}L
    OUT_RAW <- '$OUT_DIR/stability_raw${suffix}.csv'
    OUT_SUMMARY <- '$out'
    source('$ROOT/stability.R')
  " 2>&1 | tee -a "$LOG"
done

log "Done. Log: $LOG"
