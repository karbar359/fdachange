#!/usr/bin/env Rscript
# Monte-Carlo sampler for the 2-sided Hoelder-type null used by
# wilcoxon_statistic (Chapter 3 / paper 1 A1).
#
# Samples
#   T = sup_{0 <= s < t <= 1} |W(t) - W(s)| * rho(|t-s|(1-|t-s|))
# where W is a standard Brownian bridge on [0, 1] and rho is the
# symmetric weight from shared/weights.R (rho(u) = 1 / (u^gamma *
# log(c/u)^beta)).
#
# Writes the long-format CSV (x, gamma, beta) that
# shared/null_cache.R::load_wilcoxon_null reads. Existing cache
# `data-raw/cache-src/T2_gamma_beta.csv` was built with this formula (GPU-accelerated)
# at N = 10^4 grid and ~30,000 samples per (gamma, beta) cell -- weeks
# of single-node time. This pure-R port runs vectorised per lag; small
# defaults (N = 500, N_MC = 200) complete in minutes on 20 cores.
#
# Usage:
#   Rscript reproduction/generators/wilcoxon_null.R
#
# To regenerate at paper-accuracy: set N = 10000, N_MC = 30000, and
# prepare to wait a few hours per (gamma, beta) cell even in parallel.
# The limiting step is the O(N^2) per-sample pairwise-difference sweep.

# ---- Inputs (edit these) --------------------------------------------------
# Defaults are a validation-scale run (~10 min on 20 cores). For the real
# cache used by shared/null_cache.R, crank N_GRID up to 10000 and N_MC to
# 30000 -- that's hours per (gamma, beta) cell and was run on GPU.
# N_GRID = 2000, N_MC = 500 matches paper 1 Table 1 quantiles to ~0.03.
N_GRID   <- 2000L                   # lag resolution; 10000 in the real cache
N_MC     <- 500L                    # MC draws per (gamma, beta) cell
GAMMA_S  <- c(0, 0.25, 0.45)        # trimmed sweep; full sweep has 12 gammas
BETA_S_AT_GAMMA_HALF <- c(0)        # plus 1, 2 at gamma = 0.5 for the full cache
NCORES   <- max(1L, parallel::detectCores() - 2L)
SEED     <- 2L
OUT_CSV  <- NULL                    # NULL -> results/T2_gamma_beta.csv
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
source(file.path(root, "R", "legacy", "weights.R"))
source(file.path(root, "R", "legacy", "progress.R"))
suppressPackageStartupMessages(library(pvar))

# Pre-compute symmetric weights for every lag l = 1..N at the given gamma,
# beta. Vectorised: w_sym[l] = rho(l/N * (1 - l/N), gamma, beta).
holder_weight_cache <- function(N, gamma, beta) {
  ls <- seq_len(N)
  rho_gamma_beta(ls / N * (1 - ls / N), gamma, beta)
}

# One Brownian-bridge draw and corresponding sup.
#   T_max = max_l ( max_i |W[i+l] - W[i]| ) * weight[l]
holder_stat <- function(N, w_sym) {
  W <- as.numeric(pvar::rbridge(frequency = N, end = 1))
  vals <- numeric(N)
  for (l in seq_len(N)) {
    diffs <- abs(W[(l + 1L):(N + 1L)] - W[seq_len(N - l + 1L)])
    vals[l] <- max(diffs) * w_sym[l]
  }
  max(vals)
}

# Build the (gamma, beta) pairs to sweep.
cells <- list()
for (g in GAMMA_S) {
  if (isTRUE(all.equal(g, 0.5))) {
    for (b in BETA_S_AT_GAMMA_HALF) cells[[length(cells) + 1L]] <- c(g, b)
  } else {
    cells[[length(cells) + 1L]] <- c(g, 0)
  }
}

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
t0 <- Sys.time()
cat(sprintf("Wilcoxon Hoelder null: N=%d, N_MC=%d, %d cells, NCORES=%d\n",
            N_GRID, N_MC, length(cells), NCORES))

rows <- list()
for (ci in seq_along(cells)) {
  g <- cells[[ci]][1]; b <- cells[[ci]][2]
  w_sym <- holder_weight_cache(N_GRID, g, b)
  cat(sprintf("[%d/%d]  gamma=%5g beta=%g\n", ci, length(cells), g, b))
  t_cell <- Sys.time()
  xs_list <- mclapply_with_progress(
    seq_len(N_MC),
    function(i) holder_stat(N_GRID, w_sym),
    ncores = NCORES)
  xs <- unlist(xs_list)
  rows[[length(rows) + 1L]] <- data.frame(x = xs, gamma = g, beta = b)
  cat(sprintf("    q95=%.3f  cell=%s  total=%s\n\n",
              quantile(xs, 0.95),
              fmt_duration(as.numeric(difftime(Sys.time(), t_cell, units = "secs"))),
              fmt_duration(as.numeric(difftime(Sys.time(), t0,     units = "secs")))))
}

df <- do.call(rbind, rows)
if (is.null(OUT_CSV)) {
  OUT_CSV <- file.path(root, "results", "T2_gamma_beta.csv")
}
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", OUT_CSV, nrow(df)))
