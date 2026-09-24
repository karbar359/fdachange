#!/usr/bin/env Rscript
# Finite-sample (adjusted) null for Birbilas's E_{n,beta,gamma}. Sweeps
# (n, gamma, beta_id), simulates K H_0 replicates of n^{-1/2+gamma}
# E_{n,beta,gamma} per cell on FAR(1) data with no jump, and writes a
# long-format CSV that `shared/null_cache.R::load_far1_finite_null`
# reads.
#
# Implements the data-driven critical-value calibration from the
# dissertation Algorithm `alg:zero` and Equation `eq:adjusted_fit3`. The
# resulting empirical quantiles are the "adjusted-by-Bonferroni asymptotic
# critical values" Birbilas uses for n in [100, 1000].
#
# Usage:
#   Rscript reproduction/generators/far1_finite_null.R
#
# Output: code/results/E_n_beta_gamma_finite.csv. Manual
# copy-over to data-raw/cache-src/ once you have sanity-checked the quantiles:
#   cp code/results/E_n_beta_gamma_finite.csv \
#      data-raw/cache-src/E_n_beta_gamma_finite.csv

# ---- Inputs (edit these) --------------------------------------------------
.def <- function(name, value) {
  if (!exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    assign(name, value, envir = .GlobalEnv)
  }
}
.def("N_S",      c(100L, 200L, 500L))   # n sweep; add 1000 for full Birbilas grid
.def("GAMMA_S",  c(0, 0.25, 0.4, 0.45))
.def("BETA_IDS", c("beta1"))            # add beta2..beta4 for sensitivity sweep
.def("K",        500L)                  # H_0 replicates per cell
.def("L",        100L)                  # grid resolution; F4 default
.def("NCORES",   max(1L, parallel::detectCores() - 2L))
.def("SEED",     2L)
.def("OUT_CSV",  NULL)                  # NULL -> results/E_n_beta_gamma_finite.csv
# "corrected" simulates the null of far1_statistic_corrected (printed
# length^gamma denominator); the legacy cache pairs the legacy statistic.
.def("STAT_VARIANT", "legacy")          # legacy | corrected
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
source(file.path(root, "R", "legacy", "progress.R"))
suppressPackageStartupMessages(library(fdachange)); source(file.path(root, "R", "compat.R"))
source(file.path(root, "R", "legacy", "method_far1.R"))

.birb_stat <- if (identical(STAT_VARIANT, "corrected")) {
  far1_statistic_corrected
} else {
  far1_statistic
}
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
t0 <- Sys.time()
cells <- expand.grid(n = N_S, gamma = GAMMA_S, beta_id = BETA_IDS,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
cat(sprintf("Birbilas finite null: %d cells, K=%d, L=%d, NCORES=%d\n",
            nrow(cells), K, L, NCORES))

rows <- list()
for (ci in seq_len(nrow(cells))) {
  n_i <- cells$n[ci]; g_i <- cells$gamma[ci]; bid <- cells$beta_id[ci]
  beta_t <- far1_beta_registry[[bid]]
  if (is.null(beta_t)) stop("unknown beta_id: ", bid)
  cat(sprintf("[%d/%d]  n=%d gamma=%g beta=%s\n",
              ci, nrow(cells), n_i, g_i, bid))
  t_cell <- Sys.time()
  xs <- unlist(mclapply_with_progress(
    seq_len(K),
    function(k) {
      X0 <- generate_F4(n = n_i, a = 0, ec.interval = c(0, 0),
                        l = L, beta_t = beta_t)
      .birb_stat(X0, g_i, beta_t = beta_t)
    },
    ncores = NCORES))
  rows[[length(rows) + 1L]] <- data.frame(x = xs, n = n_i, gamma = g_i, beta_id = bid)
  cat(sprintf("    q95=%.3f  cell=%s  total=%s\n\n",
              quantile(xs, 0.95),
              fmt_duration(as.numeric(difftime(Sys.time(), t_cell, units = "secs"))),
              fmt_duration(as.numeric(difftime(Sys.time(), t0,     units = "secs")))))
}

df <- do.call(rbind, rows)
if (is.null(OUT_CSV)) {
  OUT_CSV <- file.path(root, "results", "E_n_beta_gamma_finite.csv")
}
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", OUT_CSV, nrow(df)))
