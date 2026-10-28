#!/usr/bin/env Rscript
# Corrected application tables: the full Table 7.1/7.2 grid recomputed with
# the 2026-09 audit corrections, on both the legacy and the fixed data
# pipeline, so every delta is attributable.
#
# Supersedes the printed application tables in four ways:
#   * corrected statistics (wilcoxon_scan_corrected, lp_scan, rkhs_scan --
#     the printed formulas; the legacy Wilcoxon scan divides by the weight it
#     should multiply by, so its printed T DECREASES in gamma);
#   * honest Ito rows: the printed Ito-Wilcoxon rows were computed with
#     gamma = 0 weights (`litgrid_hourly.R` calls get.weight.cache(n) without
#     gamma), so their gamma column was decorative; here every draw uses the
#     row's gamma, with the SAME 10 Ito paths shared across methods and
#     gammas (one seed, paths drawn once);
#   * both data pipelines: pipeline = "legacy" reproduces the print's matrix
#     (summed export overlap, end-stamp hour rotation, boundary-week
#     mislabelling); pipeline = "fixed" deduplicates, aligns hours to
#     interval starts, and orders weeks chronologically;
#   * a single window convention: every method reports I(k, m) = {k+1..m}
#     endpoints as (k/n, m/n) (the print mixed k/n and (k-1)/n across rows).
#
# Sample definitions follow the print: "with outlier" = the first 348
# calendar weeks, "no outlier" = the same span with the ISO 2021-W49 spike
# week removed. Under the fixed pipeline the samples are defined by calendar
# labels (span end 202433, outlier 202149) because the fixed matrix recovers
# two weeks the legacy one mislaid, so row INDICES no longer match.
#
# Output: code/data/litgrid_corrected_tables.csv (committed) and a LaTeX
# fragment. QUICK = TRUE trims the Ito draws for a
# sub-minute smoke run; the committed table is the QUICK = FALSE run.
#
# Usage:  Rscript applications/litgrid/corrected_tables.R

# ---- Inputs (edit these) --------------------------------------------------
if (!exists("QUICK")) QUICK <- FALSE   # pre-set QUICK <- TRUE for a smoke run
ITO_REPS     <- if (QUICK) 2L else 10L
GAMMAS_WILC  <- c(0.05, 0.25, 0.45)
GAMMAS_LP    <- c(0, 0.25, 0.45)
GAMMAS_RKHS  <- c(0, 0.25, 0.45)
ALPHA        <- 0.05
KERNEL_WILC  <- "h_0"
SEED         <- 3L
SPAN_END     <- 202433L   # last calendar week the printed 344-week trim keeps
OUTLIER_WEEK <- 202149L   # ISO 2021-W49, the removed spike week
N_WITH       <- 348L      # the print's with-outlier sample size
OUT_CSV      <- "data/litgrid_corrected_tables.csv"
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(fda); library(fda.usc) })

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) return(normalizePath(hit[1]))
  NULL
})()
root <- if (is.null(script_path)) getwd() else
  normalizePath(file.path(dirname(script_path), "..", ".."))

source(file.path(root, "R", "legacy_env.R"))
legacy_attach(root, projections = TRUE)   # vendored pre-audit pipeline
source(file.path(root, "applications", "litgrid", "data.R"))

litgrid_dir <- Sys.getenv("LITGRID_RAW_DIR", unset = "")
if (!nzchar(litgrid_dir)) {
  litgrid_dir <- normalizePath(file.path(root, "..", "data-raw", "litgrid-raw"),
                               mustWork = FALSE)
}
if (!dir.exists(litgrid_dir)) {
  stop("Litgrid raw CSVs not found at ", litgrid_dir,
       "; run Rscript applications/litgrid/download_litgrid.R or set LITGRID_RAW_DIR.")
}

# T1/T2 transforms per sec_results.tex (the framework section's chi2-only T2
# is one of the recorded errata; the results section's per-kernel pairing is
# what the tables use).
T1_uniform <- function(x) rank(x, ties.method = "average") / length(x)
T2_normal  <- function(x) qnorm(pmin(pmax(T1_uniform(x), 1e-9), 1 - 1e-9))
T2_chi2    <- function(x) qchisq(pmin(pmax(T1_uniform(x), 1e-9), 1 - 1e-9), df = 1)
RKHS_KERNEL <- list(ito = "gaussian", c_max_var = "gaussian", l2 = "min")
RKHS_BASE   <- list(ito = "",         c_max_var = "",         l2 = "rnorm.sq")
RKHS_T2     <- list(ito = T2_normal,  c_max_var = T2_normal,  l2 = T2_chi2)

load_samples <- function(pipeline) {
  if (pipeline == "legacy") {
    m <- suppressWarnings(load_litgrid_balancing_prices(
      litgrid_dir, dedup_exports = FALSE, fix_week_labels = FALSE, stamp = "to"))
    with_out <- m[seq_len(N_WITH), , drop = FALSE]
    no_out   <- m[setdiff(seq_len(nrow(m)), 204L), , drop = FALSE][seq_len(344L), , drop = FALSE]
  } else {
    m <- suppressWarnings(load_litgrid_balancing_prices(
      litgrid_dir, dedup_exports = TRUE, fix_week_labels = TRUE, stamp = "from"))
    lab <- as.integer(rownames(m))
    with_out <- m[lab <= SPAN_END, , drop = FALSE]
    no_out   <- m[lab <= SPAN_END & lab != OUTLIER_WEEK, , drop = FALSE]
  }
  list(with_outlier = with_out, no_outlier = no_out)
}

rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

week_of <- function(labels, frac, n) {
  i <- max(1L, min(n, round(frac * n)))
  labels[i]
}

for (pipeline in c("legacy", "fixed")) {
  samples <- load_samples(pipeline)
  for (sample_id in names(samples)) {
    f_matrix <- samples[[sample_id]]
    labels   <- rownames(f_matrix)
    X   <- t(f_matrix)
    n   <- ncol(X)
    Xr  <- reconstruct(X, "fs_ls", 4L, grid = "vertices")$X
    cat(sprintf("== %s / %s : n = %d ==\n", pipeline, sample_id, n))

    # The 10 Ito paths: one seed, drawn once, shared across methods and gammas.
    set.seed(SEED)
    ito_draws <- lapply(seq_len(ITO_REPS), function(r) project_ito(Xr))
    set.seed(SEED)
    proj_series <- list(
      l2        = list(project_l2(Xr)),
      c_max_var = list(project_c_max_var(Xr)),
      ito       = ito_draws)

    for (proj in names(proj_series)) {
      series <- proj_series[[proj]]
      for (method in c("wilcoxon", "lp", "rkhs")) {
        gammas <- switch(method, wilcoxon = GAMMAS_WILC,
                         lp = GAMMAS_LP, rkhs = GAMMAS_RKHS)
        for (g in gammas) {
          cv <- tryCatch(switch(method,
            wilcoxon = critical_value(
              wilcoxon_null_distribution(n, g), ALPHA, 1L),
            lp       = critical_value(lp_null_distribution(n, g, p = 2L), ALPHA, 1L),
            rkhs     = critical_value(rkhs_null_distribution(
              n, g, kernel_name = RKHS_KERNEL[[proj]],
              baseline = RKHS_BASE[[proj]], variant = "corrected"), ALPHA, 1L)),
            error = function(e) NA_real_)
          if (is.na(cv)) {
            cat(sprintf("   %-8s %-9s g=%.2f  SKIPPED (null cache not ready)\n",
                        method, proj, g))
            next
          }
          st <- lo <- hi <- numeric(length(series))
          for (b in seq_along(series)) {
            v  <- series[[b]]
            sc <- switch(method,
              wilcoxon = wilcoxon_scan_corrected(v, g, kernel = KERNEL_WILC),
              lp       = lp_scan(v, g, p = 2L, cdf = function(x) ecdf(x)(x)),
              rkhs     = rkhs_scan(RKHS_T2[[proj]](v), g,
                                   kernel = kernel_fn(RKHS_KERNEL[[proj]])))
            st[b] <- sc$stat; lo[b] <- sc$interval[1]; hi[b] <- sc$interval[2]
          }
          nd <- length(series)
          add_row(pipeline = pipeline, sample = sample_id, method = method,
                  projection = proj, gamma = g, n = n, n_draws = nd,
                  stat = mean(st), stat_sd = if (nd > 1) sd(st) else NA_real_,
                  cv = cv, reject_rate = mean(st > cv),
                  lo = mean(lo), hi = mean(hi),
                  lo_sd = if (nd > 1) sd(lo) else NA_real_,
                  hi_sd = if (nd > 1) sd(hi) else NA_real_,
                  week_lo = week_of(labels, mean(lo), n),
                  week_hi = week_of(labels, mean(hi), n))
          cat(sprintf("   %-8s %-9s g=%.2f  T=%6.3f%s  cv=%.3f  rej=%.2f  [%.3f, %.3f]\n",
                      method, proj, g, mean(st),
                      if (nd > 1) sprintf("+-%.2f", sd(st)) else "     ",
                      cv, mean(st > cv), mean(lo), mean(hi)))
        }
      }
    }
  }
}

out <- do.call(rbind, rows)
out_path <- if (QUICK) file.path(root, "results", basename(OUT_CSV)) else
  file.path(root, OUT_CSV)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_path, row.names = FALSE)
cat(sprintf("\nwrote %d rows to %s%s\n", nrow(out), out_path,
            if (QUICK) "  (QUICK run; set QUICK <- FALSE for the committed table)" else ""))
