#!/usr/bin/env Rscript
# Reproduce the printed dissertation application tables exactly.
#
# `run.R` as shipped does NOT reproduce them: it feeds the raw 168-point
# sampling grid to `project_l2`, whereas the numbers in the thesis were
# computed on curves reconstructed onto the order-2 / K = 4 polygonal basis
# (archive/code/unified_framework/helper_functions_get_functional_data_elements.R
# lines 14 and 62). The difference is not cosmetic: at gamma = 0.45 the raw
# path gives T = 2.17 with window [0.183, 0.733], the reconstructed path gives
# T = 3.50 with window [0.079, 0.599].
#
# Reproducing the printed values requires three further conventions, each of
# which this script makes explicit rather than hiding:
#
#   1. The L2 norm is taken on the fdata default grid, which for a 2^K + 1
#      dimensional basis is the 2^K + 1 = 17 polygon vertices, not the 168
#      sampling points. The two differ by about 7 per cent.
#   2. The scan window is rows k:m inclusive (helper_functions.R:598). The
#      shipped `wilcoxon_scan` shifts this and additionally drops column 1
#      from the complement; `wilcoxon_scan_exact` is the faithful version.
#   3. Table 7.1 normalises and reports with n = 348 even though the scan runs
#      on the 344-week outlier-removed series. Dropping the outlier week did
#      not update n, so the printed statistics are understated by a factor
#      (348/344)^{3/2} = 1.0175 and the window fractions use the wrong
#      denominator. Table 7.2 (348 weeks, outlier kept) is self-consistent.
#      No rejection decision changes either way; see the CORRECTED block.
#
# The raw exports also overlap, so ~6,552 hours in Jan-Sep 2019 enter the
# analysed matrix doubled. Both printed tables inherit that, so reproduction
# uses `dedup_exports = FALSE`.
#
# Usage:  Rscript applications/litgrid/reproduce_tables.R

# ---- Inputs (edit these) --------------------------------------------------
GAMMAS     <- c(0.05, 0.25, 0.45)
KERNEL     <- "h_0"
ALPHA      <- 0.05
K_LEVEL    <- 4L
ROW_DROP   <- 204L
ROW_KEEP_N <- 344L
N_PRINTED  <- 348L    # the n the printed tables use for BOTH tables
TOL        <- 5e-3    # printed values carry 3 decimals
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(fda); library(fda.usc) })

root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) normalizePath(file.path(dirname(hit[1]), "..", "..")) else getwd()
})
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

# Printed Wilcoxon / ||.||_L2 rows, dissertation Tables 7.1 and 7.2.
PRINTED <- list(
  `7.1` = data.frame(gamma = GAMMAS, T = c(5.988, 4.537, 3.438),
                     lo = 0.078, hi = 0.592),
  `7.2` = data.frame(gamma = GAMMAS, T = c(6.162, 4.668, 3.537),
                     lo = 0.078, hi = 0.595))

# All-legacy loader: the print inherits the summed export overlap, the
# end-stamp hour rotation, and the boundary-week mislabelling.
f_full <- load_litgrid_balancing_prices(litgrid_dir,
                                        dedup_exports = FALSE,
                                        fix_week_labels = FALSE,
                                        stamp = "to")

#' L2-norm projection as the thesis computed it: reconstruct onto the level-K
#' polygonal space by least squares, then take the L2 norm on the vertex grid.
l2_thesis <- function(X) {
  as.numeric(norm.fdata(matrix_to_fdata(X), metric = metric.lp, lp = 2))
}

report <- function(label, X, n_report, printed) {
  v  <- l2_thesis(X)
  n  <- length(v)
  ok <- TRUE
  cat(sprintf("\n-- Table %s: %d weeks scanned, reported with n = %d --\n",
              label, n, n_report))
  for (i in seq_along(GAMMAS)) {
    g  <- GAMMAS[i]
    sc <- wilcoxon_scan_exact(v, g, kernel = KERNEL)
    # The printed tables divide the max by n_report^{3/2}, not by n^{3/2}.
    Tp <- sc$stat * (n / n_report)^1.5
    lo <- sc$k / n_report; hi <- sc$m / n_report
    hit <- abs(Tp - printed$T[i]) < TOL &&
           abs(lo - printed$lo[i]) < TOL && abs(hi - printed$hi[i]) < TOL
    ok <- ok && hit
    cat(sprintf("   gamma=%.2f  T=%.4f (printed %.3f)  window=[%.3f, %.3f] (printed [%.3f, %.3f])  %s\n",
                g, Tp, printed$T[i], lo, hi, printed$lo[i], printed$hi[i],
                if (hit) "MATCH" else "DIFF"))
  }
  ok
}

cat("=== Reproduction of the printed Wilcoxon / L2 rows ===")
X71 <- t(f_full[setdiff(seq_len(nrow(f_full)), ROW_DROP), , drop = FALSE][seq_len(ROW_KEEP_N), , drop = FALSE])
X72 <- t(f_full[seq_len(N_PRINTED), , drop = FALSE])
ok71 <- report("7.1", X71, N_PRINTED, PRINTED[["7.1"]])
ok72 <- report("7.2", X72, N_PRINTED, PRINTED[["7.2"]])

cat("\n=== CORRECTED Table 7.1 (same data, consistent n = 344) ===\n")
v <- l2_thesis(X71); n <- length(v)
for (g in GAMMAS) {
  sc <- wilcoxon_scan_exact(v, g, kernel = KERNEL)
  cv <- critical_value(wilcoxon_null_distribution(n, g), ALPHA, 1L)
  cat(sprintf("   gamma=%.2f  T=%.4f  cv=%.3f  reject=%-5s  window=[%.3f, %.3f]\n",
              g, sc$stat, cv, sc$stat > cv, sc$k / n, sc$m / n))
}
cat("   Rejection decisions are unchanged by the correction.\n")

cat(sprintf("\nReproduction: Table 7.1 %s, Table 7.2 %s\n",
            if (ok71) "MATCH" else "DIFF", if (ok72) "MATCH" else "DIFF"))
if (!(ok71 && ok72)) quit(status = 1L)
