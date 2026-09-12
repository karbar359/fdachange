#!/usr/bin/env Rscript
# Monte-Carlo sampler for the RKHS scan-statistic null used by
# rkhs_statistic (Chapter 5 / paper 2 A2).
#
# Under H0 the sample X_1, ..., X_n is iid from the chosen baseline
# distribution (N(0,1), U(0,1), N(0,1)^2, ...). This script draws
# N_MC such samples, runs `rkhs_statistic` with the requested kernel
# and gamma, and writes the scan-statistic distribution to a wide CSV
# with columns named
#   SD.<gamma>.h.<kernel>.n<N>[.<baseline>]
# matching shared/null_cache.R::load_rkhs_null's lookup convention.
#
# The existing cache `data-raw/cache-src/SD_h_n2000_dist_t.csv` uses N = 2000,
# N_MC = 1000, three (kernel, baseline) pairs (gaussian/N(0,1),
# min/N(0,1)^2, min/U(0,1), 2_max/U(0,1)). Port uses the vectorised
# rkhs_statistic (shared with the driver); at N = 2000 one draw takes
# ~0.5s, so ~500s per cell at N_MC = 1000 in serial; ~1 minute per
# cell on 20 cores.
#
# To close the RKHS calibration gap documented in README.md (Pareto /
# Beta baselines missing), add new entries to the CELLS list.
#
# Usage: Rscript reproduction/generators/rkhs_null.R

# ---- Inputs (edit these) --------------------------------------------------
# Defaults are a fast verification run (~1 min on 20 cores). Paper 2's
# cache was N_GRID = 2000, N_MC = 1000; at those settings expect ~20 min
# per cell on 20 cores (scan is O(N^2); 12 cells -> ~4 hours).
# Constants are .def-guarded so a wrapper can pre-set them before source()
# (same pattern as driver.R); defaults unchanged.
.def <- function(name, value) {
  if (!exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    assign(name, value, envir = .GlobalEnv)
  }
}
.def("N_GRID",  150L)               # sample size (paper: 2000)
.def("N_MC",    50L)                # MC draws per cell  (paper: 1000)
.def("NCORES",  max(1L, parallel::detectCores() - 2L))
.def("SEED",    2L)
.def("OUT_CSV", NULL)               # NULL -> results/SD_h_n<N>_dist_t.csv
# "corrected" samples the null of rkhs_statistic_corrected (the printed
# formula); the legacy cache tabulates the biased legacy statistic.
.def("STAT_VARIANT", "legacy")      # legacy | corrected

# One entry per (gamma, kernel, baseline) cell to sample. Column name
# in the output CSV is built from the fields below; `baseline_tag` is
# appended as a `.<tag>` suffix (empty string -> no suffix, matching the
# existing cache's gaussian/N(0,1) column SD.<g>.h.gaussian.n2000).
CELLS <- local({
  cs <- list()
  for (g in c(0, 0.25, 0.45)) {
    cs[[length(cs) + 1L]] <- list(gamma = g, kernel = "gaussian",
                                  sampler = function(n) rnorm(n, 0, 1),
                                  baseline_tag = "")
    cs[[length(cs) + 1L]] <- list(gamma = g, kernel = "min",
                                  sampler = function(n) rnorm(n, 0, 1)^2,
                                  baseline_tag = "rnorm.sq")
    cs[[length(cs) + 1L]] <- list(gamma = g, kernel = "min",
                                  sampler = function(n) runif(n),
                                  baseline_tag = "runif")
    cs[[length(cs) + 1L]] <- list(gamma = g, kernel = "2_max",
                                  sampler = function(n) runif(n),
                                  baseline_tag = "runif")
    # The scalar S4 cell (Beta(1,9) marginal, calibrated) reads the
    # rbeta.1.9-baseline column that the legacy cache gained via
    # merge_S4_2max_beta_null.R; sample it directly here.
    cs[[length(cs) + 1L]] <- list(gamma = g, kernel = "2_max",
                                  sampler = function(n) rbeta(n, 1, 9),
                                  baseline_tag = "rbeta.1.9")
  }
  cs
})
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
source(file.path(root, "R", "legacy", "progress.R"))
source(file.path(root, "R", "legacy", "method_rkhs.R"))

.rkhs_stat <- if (identical(STAT_VARIANT, "corrected")) {
  rkhs_statistic_corrected
} else {
  rkhs_statistic
}
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
t0 <- Sys.time()
cat(sprintf("RKHS scan null (%s): N=%d, N_MC=%d, %d cells, NCORES=%d\n",
            STAT_VARIANT, N_GRID, N_MC, length(CELLS), NCORES))

columns <- list()
for (ci in seq_along(CELLS)) {
  c_ <- CELLS[[ci]]
  g <- c_$gamma; k_name <- c_$kernel; tag <- c_$baseline_tag
  kfn <- kernel_fn(k_name)
  sampler <- c_$sampler
  suffix  <- if (nzchar(tag)) paste0(".", tag) else ""
  # Column names follow the loader's convention: kernel underscores become
  # dots (load_rkhs_null looks up "2.max", as the legacy cache spells it).
  colname <- sprintf("SD.%s.h.%s.n%d%s",
                     format(g, trim = TRUE),
                     gsub("_", ".", k_name, fixed = TRUE), N_GRID, suffix)
  cat(sprintf("[%d/%d]  %s\n", ci, length(CELLS), colname))
  t_cell <- Sys.time()
  xs_list <- mclapply_with_progress(
    seq_len(N_MC),
    function(i) .rkhs_stat(sampler(N_GRID), gamma = g, kernel = kfn),
    ncores = NCORES)
  xs <- unlist(xs_list)
  columns[[colname]] <- xs
  cat(sprintf("    q95=%.4f  cell=%s  total=%s\n\n",
              quantile(xs, 0.95),
              fmt_duration(as.numeric(difftime(Sys.time(), t_cell, units = "secs"))),
              fmt_duration(as.numeric(difftime(Sys.time(), t0,     units = "secs")))))
}

df <- as.data.frame(columns, check.names = FALSE)
if (is.null(OUT_CSV)) {
  OUT_CSV <- file.path(root, "results",
                       sprintf("SD_h_n%d_dist_t.csv", N_GRID))
}
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("wrote %s (%d rows x %d cols)\n", OUT_CSV, nrow(df), ncol(df)))
