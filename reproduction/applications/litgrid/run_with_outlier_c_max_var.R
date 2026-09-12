#!/usr/bin/env Rscript
# Antisymmetric scan with C-max-var projection on the FULL (outlier-included)
# Litgrid balancing-price dataset, evaluated at gamma in {0.05, 0.25, 0.45}.
# Mirrors run.R but skips the ROW_DROP step so the outlier week is retained.
# Prints a TSV-ish summary that can be pasted into tab:litgrid_asym.

KERNEL     <- "h_0"
ALPHA      <- 0.05
ROW_KEEP_N <- 344L            # keep the first 344 valid weeks
SEED       <- 3L

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
f_matrix <- load_litgrid_balancing_prices(litgrid_dir)
cat(sprintf("Raw weekly matrix: %d weeks x %d hours\n",
            nrow(f_matrix), ncol(f_matrix)))

# Keep first ROW_KEEP_N rows AS-IS (outlier week 204 retained).
f_matrix <- f_matrix[seq_len(ROW_KEEP_N), , drop = FALSE]
cat(sprintf("After trim (no row drop, keep %d): %d weeks x %d hours\n",
            ROW_KEEP_N, nrow(f_matrix), ncol(f_matrix)))

# Orient rows = grid (168 hours), cols = observations (n weeks).
X <- t(f_matrix)

set.seed(SEED)
project <- projection_fn("c_max_var")
X_proj  <- project(X)
cat(sprintf("Projected to numeric(%d) via c_max_var (range=[%.3g, %.3g], sd=%.3g)\n",
            length(X_proj), min(X_proj), max(X_proj), sd(X_proj)))

cat("\ngamma\tT\tk_idx\tm_idx\tlow\thigh\tcv95\treject\n")
for (gamma in c(0.05, 0.25, 0.45)) {
  scan <- wilcoxon_scan(X_proj, gamma, kernel = KERNEL)
  null <- wilcoxon_null_distribution(ncol(X), gamma)
  cv   <- critical_value(null, ALPHA, m.d = 1L)
  cat(sprintf("%.2f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%s\n",
              gamma, scan$stat, scan$k, scan$m,
              scan$interval[1], scan$interval[2], cv,
              scan$stat > cv))
}
