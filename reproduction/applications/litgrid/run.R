#!/usr/bin/env Rscript
# Litgrid practical application (paper 1 §6, dissertation Ch6).
#
# Loads Lithuanian balancing-energy-price data from <repo>/.raw/litgrid/,
# projects each week-long curve down to a scalar via the chosen functional
# projection, and runs the Wilcoxon scan test on the resulting weekly
# series to detect an epidemic change segment. Prints the detected
# (start_week, end_week) interval with the associated test statistic and
# 5% critical value from the Hoelder-type null cache.
#
# Default configuration matches the dissertation (L2 projection, h_0
# Wilcoxon kernel, gamma = 0.45).
#
# Date erratum (2026-09-05): the thesis and older comments here tied the
# removed outlier week to the 2021-11-20 grid disturbance (ISO 2021-W46).
# Verified from the raw timestamps, the removed row 204 is ISO 2021-W49,
# 6-12 December 2021, whose balancing prices peak at 1199.9 EUR/MWh on
# 2021-12-08 09:00; week 2021-W46 peaks at 309.5 and is unremarkable.
#
# Usage:  Rscript applications/litgrid/run.R
# Env:    LITGRID_RAW_DIR overrides the raw-data lookup (default:
#         <repo>/data-raw/litgrid-raw, fetched by download_litgrid.R).

# ---- Inputs (edit these) --------------------------------------------------
GAMMA      <- 0.45                  # Hoelder-weight; paper 1 reports .45
PROJECTION <- "l2"                  # l2 | ito | wiener | c_max_var | pc1 | pc_d
KERNEL     <- "h_0"                 # h_0 (norm) -- paper 1's choice
ALPHA      <- 0.05
# Curve reconstruction, applied before the projection. The dissertation
# smooths each 168-point week onto the level-K polygonal (Faber-Schauder /
# order-2 B-spline) space; omitting this step does NOT reproduce the printed
# tables. See the reconstruct() documentation in the fdachange package.
RECONSTRUCTION <- "fs_ls"           # none | fs_interp | fs_ls | bspline_pen
RECON_K        <- 4L                # dissertation Sec. 7.2 uses K = 4
RECON_EVAL     <- "vertices"        # vertices (2^K+1 points) | input
DEDUP_EXPORTS  <- TRUE              # FALSE reproduces the printed tables
FIX_WEEK_LABELS <- TRUE             # FALSE keeps the boundary-week mislabel
STAMP          <- "from"            # "to" keeps the one-hour curve rotation
# Scan family: "corrected" is the printed statistic (multiply-weight, true
# window length, Kuiper gamma = 0 null); "exact" fixes only the summed-area
# off-by-one and still carries the inverted weight; "legacy" is the scan that
# fed the print.
SCAN           <- "corrected"       # corrected | exact | legacy
# Sample trim: drop the outlier week (ISO 2021-W49, the 6-12 Dec 2021
# price-spike week; a genuine market event, not a data artefact) and keep the
# calendar span the printed 344-week trim covered. Selection is BY LABEL, not
# by row index: under the fixed pipeline the outlier sits at row 206, not the
# legacy 204, because two mislaid weeks are recovered.
DROP_WEEK  <- 202149L               # ISO label of the outlier week
SPAN_END   <- 202433L               # last week the printed trim kept
SEED       <- 3L
# ---------------------------------------------------------------------------

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
  normalizePath(file.path(dirname(script_path), "..", ".."))
}

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
litgrid_dir <- litgrid_dir

cat(sprintf("Loading Litgrid balancing prices from %s\n", litgrid_dir))
f_matrix <- load_litgrid_balancing_prices(litgrid_dir,
                                          dedup_exports = DEDUP_EXPORTS,
                                          fix_week_labels = FIX_WEEK_LABELS,
                                          stamp = STAMP)
cat(sprintf("Raw weekly matrix: %d weeks x %d hours\n",
            nrow(f_matrix), ncol(f_matrix)))

# Label-based trim (robust to the pipeline recovering mislaid weeks).
lab <- as.integer(rownames(f_matrix))
f_matrix <- f_matrix[lab <= SPAN_END & lab != DROP_WEEK, , drop = FALSE]
cat(sprintf("After trim (drop %d, span <= %d): %d weeks x %d hours\n",
            DROP_WEEK, SPAN_END, nrow(f_matrix), ncol(f_matrix)))

# Orient as (l+1 = 168 hours) rows x n = 344 weeks cols -- matches
# code's functional-scenario contract (rows = grid,
# cols = observations).
X <- t(f_matrix)

set.seed(SEED)
rec <- reconstruct(X, RECONSTRUCTION, RECON_K, grid = RECON_EVAL)
cat(sprintf("Reconstruction: %s (grid %d points)\n",
            recon_label(RECONSTRUCTION, RECON_K, RECON_EVAL), nrow(rec$X)))
project <- projection_fn(PROJECTION)
X_proj  <- project(rec$X)
cat(sprintf("Projected to numeric(%d) via %s (range=[%.3g, %.3g], sd=%.3g)\n",
            length(X_proj), PROJECTION, min(X_proj), max(X_proj), sd(X_proj)))

scan <- switch(SCAN,
  corrected = wilcoxon_scan_corrected(X_proj, GAMMA, kernel = KERNEL),
  exact     = wilcoxon_scan_exact(X_proj, GAMMA, kernel = KERNEL),
  legacy    = wilcoxon_scan(X_proj, GAMMA, kernel = KERNEL),
  stop("SCAN must be corrected | exact | legacy"))
null <- wilcoxon_null_distribution(ncol(X), GAMMA)
cv   <- critical_value(null, ALPHA, m.d = 1L)

cat("\n==== result ====\n")
cat(sprintf("n=%d  gamma=%.2f  projection=%s  kernel=%s\n",
            ncol(X), GAMMA, PROJECTION, KERNEL))
cat(sprintf("stat  = %.4f\n", scan$stat))
cat(sprintf("cv(95%%) = %.4f\n", cv))
cat(sprintf("reject H0 (no change) at alpha=%.2f: %s\n",
            ALPHA, scan$stat > cv))
cat(sprintf("detected window (column indices)  : [k=%d, m=%d]\n",
            scan$k, scan$m))
cat(sprintf("detected window (fraction of n)   : [%.3f, %.3f]\n",
            scan$interval[1], scan$interval[2]))

# Week labels (row names preserved in prep.data.into.matrix as year_week).
week_labels <- rownames(f_matrix)
if (length(week_labels) >= scan$m) {
  cat(sprintf("detected window (ISO year-week)   : %s .. %s\n",
              week_labels[scan$k], week_labels[scan$m]))
}
