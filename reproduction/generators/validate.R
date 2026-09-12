#!/usr/bin/env Rscript
# Compare the CSVs under results/ (produced by the
# three null generators) against the canonical caches in data-raw/cache-src/.
#
# Read-only: picks up whatever the generators last wrote to results/,
# matches columns to the cache, and prints a quantile / KS summary. If
# you want a stricter match, first rerun the relevant generator at
# settings closer to the cache's (see each generator's header comment
# for paper-accuracy parameters), then rerun this script.
#
# Column matching:
#   Wilcoxon : results/T2_gamma_beta.csv      vs  cache-src/T2_gamma_beta.csv
#              (match on (gamma, beta) rows; cache was built at N=10^4;
#               my generator runs at a smaller N so a finite-N shift in
#               quantiles is expected, not a bug).
#   Lp       : results/dist_p2_ws_n*_t.csv    vs  cache-src/dist_p2_ws_n200_t.csv
#              (match columns by name; N-suffix must match for a fair
#               comparison -- warns if not).
#   RKHS     : results/SD_h_n*_dist_t.csv     vs  cache-src/SD_h_n2000_dist_t.csv
#              (same column-name match; N likewise must match).
#
# Usage: Rscript reproduction/generators/validate.R
# REPRO_CACHE_DIR env var overrides the cache lookup.

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- sub(key, "", args[grepl(key, args)])
  if (length(hit)) return(normalizePath(hit[1]))
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) return(normalizePath(ofile))
  NULL
})()
root <- if (is.null(script_path)) {
  getwd()
} else {
  normalizePath(file.path(dirname(script_path), ".."))
}

raw_dir <- Sys.getenv("REPRO_CACHE_DIR", unset = NA)
if (is.na(raw_dir) || !nzchar(raw_dir)) {
  cur <- root; raw_dir <- NULL
  for (i in 1:6) {
    cand <- file.path(cur, "data-raw", "cache-src")
    if (dir.exists(cand)) { raw_dir <- cand; break }
    cur <- dirname(cur)
  }
  if (is.null(raw_dir)) stop("Could not find data-raw/cache-src/; set REPRO_CACHE_DIR.")
}
results_dir <- file.path(root, "results")

# Published quantiles from paper 1 (A1_LMJ.tex Table 1 `tab:T` for T_gamma
# with beta = 0, and Table 2 `tab:T-beta` for gamma = 0.5). Paper generated
# these at grid size 10^4 with 30000 MC runs, so they serve as the
# ground-truth reference both our generator and the cache must match.
# Column order: quantile levels 50%, 80%, 90%, 95%, 99%.
PAPER_T_GAMMA <- list(
  "0"     = c(1.213, 1.460, 1.612, 1.741, 2.012),
  "0.05"  = c(1.314, 1.573, 1.732, 1.862, 2.143),
  "0.1"   = c(1.423, 1.708, 1.876, 2.016, 2.306),
  "0.15"  = c(1.544, 1.839, 2.017, 2.172, 2.485),
  "0.2"   = c(1.684, 1.995, 2.175, 2.344, 2.677),
  "0.25"  = c(1.846, 2.172, 2.357, 2.527, 2.862),
  "0.3"   = c(2.042, 2.372, 2.572, 2.748, 3.122),
  "0.35"  = c(2.294, 2.621, 2.825, 3.017, 3.387),
  "0.4"   = c(2.654, 2.961, 3.150, 3.330, 3.695),
  "0.45"  = c(3.268, 3.529, 3.697, 3.852, 4.216))

summarise_pair <- function(mine, cache, label, matched_N = TRUE) {
  if (length(mine) == 0L || length(cache) == 0L) {
    cat(sprintf("  %-40s  SKIP (len mine=%d, cache=%d)\n",
                label, length(mine), length(cache)))
    return(invisible())
  }
  qs <- c(0.05, 0.5, 0.95, 0.99)
  qm <- quantile(mine,  qs); qc <- quantile(cache, qs)
  ks <- suppressWarnings(ks.test(mine, cache))
  # MC standard-error rough proxy: sd / sqrt(length). Gives a sense of how
  # tight each quantile estimate is.
  verdict <- if (!matched_N) {
    "N-MISMATCH"
  } else if (ks$p.value > 0.05) {
    "match"
  } else if (ks$p.value > 0.01) {
    "borderline"
  } else {
    "DIFF"
  }
  cat(sprintf(
    "  %-40s  n=%4d/%-5d  q05=%5.3f/%5.3f  q50=%5.3f/%5.3f  q95=%5.3f/%5.3f  q99=%5.3f/%5.3f  mean=%5.3f/%5.3f  sd=%5.3f/%5.3f  KSp=%.3f  %s\n",
    label, length(mine), length(cache),
    qm[1], qc[1], qm[2], qc[2], qm[3], qc[3], qm[4], qc[4],
    mean(mine), mean(cache), sd(mine), sd(cache),
    ks$p.value, verdict))
}

# ==== Wilcoxon ============================================================
wilc_path <- file.path(results_dir, "T2_gamma_beta.csv")
wilc_cache_path <- file.path(raw_dir, "T2_gamma_beta.csv")
if (file.exists(wilc_path) && file.exists(wilc_cache_path)) {
  cat("==== Wilcoxon Hoelder null ====\n")
  cat(sprintf("  mine : %s\n  cache: %s\n",
              wilc_path, wilc_cache_path))
  mine  <- read.csv(wilc_path)
  cache <- read.csv(wilc_cache_path)
  cat(sprintf(
    "  NOTE: cache was built at N=10^4 on GPU; generator default is N<<10^4,\n"))
  cat(sprintf(
    "        so finite-N bias is expected (q95 systematically lower than cache).\n\n"))
  # Match on (gamma, beta).
  cells <- unique(mine[, c("gamma", "beta")])
  for (i in seq_len(nrow(cells))) {
    g <- cells$gamma[i]; b <- cells$beta[i]
    xs_mine  <- mine$x[mine$gamma == g & mine$beta == b]
    xs_cache <- cache$x[cache$gamma == g & cache$beta == b]
    summarise_pair(xs_mine, xs_cache,
                   sprintf("T2 gamma=%s beta=%s", g, b),
                   matched_N = FALSE)   # N never matches for Wilcoxon
    # Also compare each against paper 1's Table 1 quantiles if available.
    key <- format(g, trim = TRUE)
    if (b == 0 && !is.null(PAPER_T_GAMMA[[key]])) {
      pq <- PAPER_T_GAMMA[[key]]
      qs <- c(0.5, 0.8, 0.9, 0.95, 0.99)
      qm <- quantile(xs_mine,  qs)
      qc <- quantile(xs_cache, qs)
      cat(sprintf("    paper1 Tab1: q50=%.3f q80=%.3f q90=%.3f q95=%.3f q99=%.3f\n",
                  pq[1], pq[2], pq[3], pq[4], pq[5]))
      cat(sprintf("    mine  - p1:  %+.3f %+.3f %+.3f %+.3f %+.3f\n",
                  qm[1] - pq[1], qm[2] - pq[2], qm[3] - pq[3],
                  qm[4] - pq[4], qm[5] - pq[5]))
      cat(sprintf("    cache - p1:  %+.3f %+.3f %+.3f %+.3f %+.3f\n",
                  qc[1] - pq[1], qc[2] - pq[2], qc[3] - pq[3],
                  qc[4] - pq[4], qc[5] - pq[5]))
    }
  }
  cat("\n")
} else {
  cat("[skip] Wilcoxon: results/T2_gamma_beta.csv not found\n",
      "       (run generators/wilcoxon_null.R to produce it)\n\n")
}

# ==== Lp ==================================================================
lp_cache_path <- file.path(raw_dir, "dist_p2_ws_n200_t.csv")
lp_results <- list.files(results_dir,
                         pattern = "^dist_p2_ws_n\\d+_t\\.csv$",
                         full.names = TRUE)
if (length(lp_results) && file.exists(lp_cache_path)) {
  cache <- read.csv(lp_cache_path, check.names = FALSE)
  for (lp_path in lp_results) {
    n_match <- regmatches(basename(lp_path),
                          regexpr("n\\d+", basename(lp_path)))
    mine_N <- as.integer(sub("n", "", n_match))
    cat(sprintf("==== Lp Wiener-sheet null (mine N=%d, cache N=200) ====\n",
                mine_N))
    cat(sprintf("  mine : %s\n  cache: %s\n\n", lp_path, lp_cache_path))
    mine <- read.csv(lp_path, check.names = FALSE)
    # My columns name N=mine_N, cache names N=200. Match by (gamma, p).
    for (col_c in names(cache)) {
      stem <- sub("\\.n\\d+\\.", sprintf(".n%d.", mine_N), col_c)
      if (stem %in% names(mine)) {
        summarise_pair(mine[[stem]], cache[[col_c]], col_c,
                       matched_N = (mine_N == 200L))
      } else {
        cat(sprintf("  %-40s  (mine missing column %s)\n", col_c, stem))
      }
    }
    cat("\n")
  }
} else {
  cat("[skip] Lp: no results/dist_p2_ws_n*_t.csv found\n\n")
}

# ==== RKHS ================================================================
rkhs_cache_path <- file.path(raw_dir, "SD_h_n2000_dist_t.csv")
rkhs_results <- list.files(results_dir,
                           pattern = "^SD_h_n\\d+_dist_t\\.csv$",
                           full.names = TRUE)
if (length(rkhs_results) && file.exists(rkhs_cache_path)) {
  cache <- read.csv(rkhs_cache_path, check.names = FALSE)
  for (rkhs_path in rkhs_results) {
    n_match <- regmatches(basename(rkhs_path),
                          regexpr("n\\d+", basename(rkhs_path)))
    mine_N <- as.integer(sub("n", "", n_match))
    cat(sprintf("==== RKHS scan null (mine N=%d, cache N=2000) ====\n",
                mine_N))
    cat(sprintf("  mine : %s\n  cache: %s\n\n", rkhs_path, rkhs_cache_path))
    mine <- read.csv(rkhs_path, check.names = FALSE)
    for (col_c in names(cache)) {
      stem <- sub("\\.n2000\\b", sprintf(".n%d", mine_N), col_c)
      if (stem %in% names(mine)) {
        summarise_pair(mine[[stem]], cache[[col_c]], col_c,
                       matched_N = (mine_N == 2000L))
      } else {
        cat(sprintf("  %-40s  (mine missing column %s)\n", col_c, stem))
      }
    }
    cat("\n")
  }
} else {
  cat("[skip] RKHS: no results/SD_h_n*_dist_t.csv found\n\n")
}
