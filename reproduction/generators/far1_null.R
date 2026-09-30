#!/usr/bin/env Rscript
# Monte-Carlo sampler for the asymptotic null E_{q,gamma} of Birbilas's
# epidemic-change statistic E_{n,beta,gamma} (`methods/birbilas/run.R`).
#
# Discretisation matches the reference R script
# the original (GPU-parallelised) table-generation script:
# build an N x N matrix whose rows are independent discretized standard
# Brownian-motion paths (cumsum of N iid N(0, 1)). For every interval
# [m, k] with k > m take max over t of |colSums(x1[m:k, ])| / (k - m)^gamma.
# Also include the single-row max(|x1[k, ]|) terms (boundary contribution).
# Multiply the overall max by N^{-1 + gamma}.
#
# Defaults are sub-minute verification scale; bump to N = 1200, N_MC = 1000
# to match Birbilas's full-fidelity setting (a few minutes per gamma on
# 20 cores, one cache per gamma cell).
#
# Usage:
#   Rscript reproduction/generators/far1_null.R

# ---- Inputs (edit these) --------------------------------------------------
N        <- 200L                     # grid resolution; 1200 in the real cache
N_MC     <- 50L                      # MC draws per gamma cell; 1000 paper-scale
GAMMA_S  <- c(0, 0.25, 0.4, 0.45)    # 0.4 is Birbilas's default for E_{n,beta,gamma}
NCORES   <- max(1L, parallel::detectCores() - 2L)
SEED     <- 2L
OUT_CSV  <- NULL                     # NULL -> results/E_q_gamma.csv
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

# Single E_{q, gamma} draw at grid resolution N.
e_q_gamma_stat <- function(N, gamma) {
  # x1: rows are BM paths. Generate as cumsum-of-rnorm per row.
  x1 <- matrix(rnorm(N * N), nrow = N, ncol = N)
  x1 <- t(apply(x1, 1L, cumsum))                # N x N, each row a BM path
  # Column prefix sums (over rows): CS[i + 1, ] = sum of x1[1..i, ] per col.
  CS <- rbind(0, apply(x1, 2L, cumsum))         # (N + 1) x N

  E_max <- 0
  # Interval part: pairs (m, k) with 1 <= m <= k <= N, sum rows m..k
  # (k - m + 1 rows), divide by (k - m)^gamma. d = k - m in 1..N-1.
  for (d in seq_len(N - 1L)) {
    m_idx <- seq.int(1L, N - d)                  # m in 1..N - d
    k_idx <- m_idx + d                            # k = m + d
    # CS[k + 1, ] - CS[m, ] = sum of x1[m..k, ] = (k - m + 1) rows.
    block <- abs(CS[k_idx + 1L, , drop = FALSE] - CS[m_idx, , drop = FALSE])
    val_d <- max(block) / d^gamma
    if (val_d > E_max) E_max <- val_d
  }
  # Single-row terms: max_t |x1[k, t]| over all (k, t).
  row_max <- max(abs(x1))
  if (row_max > E_max) E_max <- row_max
  E_max * N^(-1 + gamma)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
t0 <- Sys.time()
cat(sprintf("Birbilas E_{q,gamma} null: N=%d, N_MC=%d, %d gammas, NCORES=%d\n",
            N, N_MC, length(GAMMA_S), NCORES))

rows <- list()
for (gi in seq_along(GAMMA_S)) {
  g <- GAMMA_S[gi]
  cat(sprintf("[%d/%d]  gamma=%g\n", gi, length(GAMMA_S), g))
  t_cell <- Sys.time()
  xs <- unlist(mclapply_with_progress(
    seq_len(N_MC),
    function(i) e_q_gamma_stat(N, g),
    ncores = NCORES))
  rows[[length(rows) + 1L]] <- data.frame(x = xs, gamma = g)
  cat(sprintf("    q95=%.3f  cell=%s  total=%s\n\n",
              quantile(xs, 0.95),
              fmt_duration(as.numeric(difftime(Sys.time(), t_cell, units = "secs"))),
              fmt_duration(as.numeric(difftime(Sys.time(), t0,     units = "secs")))))
}

df <- do.call(rbind, rows)
if (is.null(OUT_CSV)) {
  OUT_CSV <- file.path(root, "results", "E_q_gamma.csv")
}
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", OUT_CSV, nrow(df)))
