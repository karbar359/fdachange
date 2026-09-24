#!/usr/bin/env Rscript
# Finite-sample null for the RKHS scan statistic on F4 (FAR(1) + Wiener
# innovations) data. Sweeps a hand-picked list of (projection, kernel)
# pairs that align with the kernel's natural domain:
#
#   (c_max_var, gaussian)  -- real-valued projection, generic kernel
#   (ito,       gaussian)  -- real-valued projection, generic kernel
#   (l2,        min)       -- non-negative projection, BM-covariance kernel
#                              (matches Ch.5 K3 design where `min` is paired
#                              with chi-squared data).
#
# Sweep over n in N_S, gamma in GAMMA_S, the three pairs above. K H_0
# replicates per cell on FAR(1) data with no jump. Output is long-format
# CSV that `shared/null_cache.R::load_rkhs_F4_null` reads.
#
# Usage:
#   Rscript reproduction/generators/rkhs_F4_null.R
#
# Output: code/results/SD_h_F4_finite.csv. Manual copy-over
# to data-raw/cache-src/ once you have sanity-checked the quantiles:
#   cp code/results/SD_h_F4_finite.csv \
#      data-raw/cache-src/SD_h_F4_finite.csv

# ---- Inputs (edit these) --------------------------------------------------
.def <- function(name, value) {
  if (!exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    assign(name, value, envir = .GlobalEnv)
  }
}
.def("N_S",     c(200L))                 # n grid; matches the F4 paper-grade run
.def("GAMMA_S", c(0, 0.25, 0.45))
.def("CELLS", list(
  c(projection = "c_max_var", kernel = "gaussian"),
  c(projection = "ito",       kernel = "gaussian"),
  c(projection = "l2",        kernel = "min")))
.def("K",       1000L)                   # H_0 replicates per cell
.def("L",       100L)                    # grid resolution; F4 default
.def("NCORES",  max(1L, parallel::detectCores() - 2L))
.def("SEED",    2L)
.def("OUT_CSV", NULL)                    # NULL -> results/SD_h_F4_finite.csv
# "corrected" simulates the null of rkhs_statistic_corrected (printed formula).
.def("STAT_VARIANT", "legacy")           # legacy | corrected
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
  normalizePath(file.path(dirname(script_path), ".."))
}
source(file.path(root, "R", "legacy", "integral_image.R"))
source(file.path(root, "R", "legacy", "weights.R"))
source(file.path(root, "R", "legacy", "kernels.R"))
source(file.path(root, "R", "legacy", "projections.R"))
source(file.path(root, "R", "legacy", "progress.R"))
suppressPackageStartupMessages(library(fdachange)); source(file.path(root, "R", "compat.R"))
source(file.path(root, "R", "legacy", "method_rkhs.R"))
source(file.path(root, "R", "legacy", "method_far1.R"))         # for beta registry

beta_t <- far1_beta_registry[["beta1"]]                      # F4 default

.rkhs_stat <- if (identical(STAT_VARIANT, "corrected")) {
  rkhs_statistic_corrected
} else {
  rkhs_statistic
}
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
t0 <- Sys.time()
plan <- expand.grid(n = N_S, gamma = GAMMA_S, cell_idx = seq_along(CELLS),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
cat(sprintf("RKHS F4 finite null: %d cells, K=%d, L=%d, NCORES=%d\n",
            nrow(plan), K, L, NCORES))

rows <- list()
for (i in seq_len(nrow(plan))) {
  n_i  <- plan$n[i]; g_i <- plan$gamma[i]; ci <- plan$cell_idx[i]
  pname  <- CELLS[[ci]][["projection"]]
  kname  <- CELLS[[ci]][["kernel"]]
  proj   <- projection_fn(pname)
  K_fn   <- kernel_fn(kname)
  cat(sprintf("[%d/%d]  n=%d gamma=%g projection=%s kernel=%s\n",
              i, nrow(plan), n_i, g_i, pname, kname))
  t_cell <- Sys.time()
  xs <- unlist(mclapply_with_progress(
    seq_len(K),
    function(k) {
      X  <- generate_F4(n = n_i, a = 0, ec.interval = c(0, 0),
                        l = L, beta_t = beta_t)
      Xp <- proj(X)
      .rkhs_stat(Xp, gamma = g_i, kernel = K_fn)
    },
    ncores = NCORES))
  rows[[length(rows) + 1L]] <- data.frame(
    x = xs, n = n_i, gamma = g_i, projection = pname, kernel = kname)
  cat(sprintf("    q95=%.3f  cell=%s  total=%s\n\n",
              quantile(xs, 0.95),
              fmt_duration(as.numeric(difftime(Sys.time(), t_cell, units = "secs"))),
              fmt_duration(as.numeric(difftime(Sys.time(), t0,     units = "secs")))))
}

df <- do.call(rbind, rows)
if (is.null(OUT_CSV)) {
  OUT_CSV <- file.path(root, "results", "SD_h_F4_finite.csv")
}
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", OUT_CSV, nrow(df)))
