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
# `data-raw/cache-src/T2_gamma_beta.csv` holds the gamma > 0 cells from the
# original GPU run (N = 10^4 grid, ~30,000 samples per (gamma, beta) cell)
# and the gamma = 0 cell regenerated with this port on 2026-10-28 at the
# same scale (tag `port_wilcoxon_null_N10000_MC30000_independent`). This
# pure-R port runs vectorised per lag: one N = 10^4 draw costs about 0.3 s,
# so a 30,000-draw cell takes about 10 minutes on 16 cores.
#
# Usage:
#   Rscript reproduction/generators/wilcoxon_null.R
#
# To regenerate a cell at paper scale, pre-set the constants before
# source(): N_GRID <- 10000L; N_MC <- 30000L; GAMMA_S <- 0; NCORES <- 16L.
# The limiting step is the O(N^2) per-sample pairwise-difference sweep.

# ---- Inputs (edit these) --------------------------------------------------
# Constants are .def-guarded so a wrapper can pre-set them before source()
# (same pattern as the other generators). Defaults are a validation-scale
# run: N_GRID = 2000, N_MC = 500 matches paper 1 Table 1 quantiles to ~0.03.
.def <- function(name, value) {
  if (!exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    assign(name, value, envir = .GlobalEnv)
  }
}
.def("N_GRID",   2000L)                 # lag resolution; 10000 in the real cache
.def("N_MC",     500L)                  # MC draws per (gamma, beta) cell
.def("GAMMA_S",  c(0, 0.25, 0.45))      # trimmed sweep; full sweep has 12 gammas
.def("BETA_S_AT_GAMMA_HALF", c(0))      # plus 1, 2 at gamma = 0.5 for the full cache
.def("NCORES",   max(1L, parallel::detectCores() - 2L))
.def("SEED",     2L)
.def("OUT_CSV",  NULL)                  # NULL -> results/T2_gamma_beta.csv
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
