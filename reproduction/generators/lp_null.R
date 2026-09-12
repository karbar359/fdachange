#!/usr/bin/env Rscript
# Monte-Carlo sampler for the Lp Wiener-sheet null used by lp_statistic
# (Chapter 4 / paper 2 A2).
#
# Samples
#   T = max_{t != s} ( integral_0^1 |delta(u, t, s)|^p du )^(1/p)
#        * rho(|t-s|(1-|t-s|))
# with
#   delta(u, t, s) = W(u, t) - W(u, s) - u*W(1, t) - (t-s)*W(u, 1)
#                  + u*W(1, s) + (t-s)*u*W(1, 1)
# the centered Wiener-sheet increment on (t, s) in [0, 1]^2, and rho the
# symmetric weight from shared/weights.R.
#
# Writes the wide-format CSV expected by shared/null_cache.R::load_lp_null:
# columns `dist.p2.ws.g<gamma>.n<N>.p<p>` (one per (gamma, p) cell). The
# cached `data-raw/cache-src/dist_p2_ws_n200_t.csv` uses N = 200, N_MC = 300, gamma in
# {0, 0.25, 0.45}, p in {2, 4, 10}. Port uses the same defaults.
#
# Usage: Rscript reproduction/generators/lp_null.R

# ---- Inputs (edit these) --------------------------------------------------
# Defaults are a fast verification run (a few minutes on 20 cores).
# Paper 2's cache uses N_GRID = 200, N_MC = 300 and takes ~4 minutes here.
N_GRID   <- 100L                    # Wiener-sheet grid size (paper 2: 200)
N_MC     <- 100L                    # MC draws per (gamma, p) cell
GAMMA_S  <- c(0, 0.25, 0.45)
P_S      <- c(2L, 4L, 10L)
NCORES   <- max(1L, parallel::detectCores() - 2L)
SEED     <- 2L
OUT_CSV  <- NULL                    # NULL -> results/dist_p2_ws_n<N>_t.csv
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
source(file.path(root, "R", "legacy", "integral_image.R"))
source(file.path(root, "R", "legacy", "progress.R"))

# Compute the Wiener-sheet sample W(i, j) for (i, j) in {1..N}^2, where
# W(i, j) = sum of N(0, 1/N^2) increments up to (i, j). Built via the
# integral image of a random NxN matrix of those iid increments.
wiener_sheet <- function(N) {
  Z <- matrix(rnorm(N * N, 0, 1), N, N) / N
  NumMatrix(Z)$sat[2:(N + 1L), 2:(N + 1L)]
}

# Weight cache indexed by lag l in 1..(N-1). Matches the legacy
# `get.weight.cached` (non-symmetric; same as methods/lp/run.R):
#   w[l] = rho(l / N, gamma, beta)
lp_weight_cache <- function(N, gamma, beta = 0) {
  ls <- seq_len(N - 1L)
  rho_gamma_beta(ls / N, gamma, beta)
}

# One (N, gamma, p) statistic draw.
#   For every ordered pair (k, l) with k != l in {1..N}, compute
#   I_{k,l} = (mean over u of |delta|^p)^(1/p) and weight by w[|k-l|].
#   Return max over pairs.
lp_sheet_stat <- function(N, gamma, p, w) {
  W <- wiener_sheet(N)
  W_end_row <- W[N, ]
  W_end_col <- W[, N]
  W_end_end <- W[N, N]
  us <- seq_len(N) / N

  # For efficiency, enumerate only pairs with k < l (symmetric in the
  # formula up to sign, and the |.|^p is sign-insensitive). Halves the
  # work compared to all ordered pairs.
  kl <- which(upper.tri(matrix(FALSE, N, N), diag = FALSE), arr.ind = TRUE)
  k_v <- kl[, 1L]; l_v <- kl[, 2L]
  lag_v <- l_v - k_v

  best <- 0
  for (p_idx in seq_along(k_v)) {
    k <- k_v[p_idx]; l <- l_v[p_idx]
    # delta_u for u = 1..N
    delta <- (W[, l] - W[, k]) - (us * W_end_row[l]) -
             ((l - k) / N) * W_end_col + (us * W_end_row[k]) +
             ((l - k) / N) * us * W_end_end
    int_u <- sum(abs(delta)^p) / N
    val   <- int_u^(1 / p) * w[lag_v[p_idx]]
    if (val > best) best <- val
  }
  best
}

# Flat (gamma, p) task list.
cells <- list()
for (g in GAMMA_S) for (p in P_S) cells[[length(cells) + 1L]] <- list(gamma = g, p = p)

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
t0 <- Sys.time()
cat(sprintf("Lp Wiener-sheet null: N=%d, N_MC=%d, %d cells, NCORES=%d\n",
            N_GRID, N_MC, length(cells), NCORES))

columns <- list()
for (ci in seq_along(cells)) {
  cell <- cells[[ci]]
  g <- cell$gamma; p <- cell$p
  w <- lp_weight_cache(N_GRID, g, beta = 0)
  cat(sprintf("[%d/%d]  gamma=%.2f p=%2d\n", ci, length(cells), g, p))
  t_cell <- Sys.time()
  xs_list <- mclapply_with_progress(
    seq_len(N_MC),
    function(i) lp_sheet_stat(N_GRID, g, p, w),
    ncores = NCORES)
  xs <- unlist(xs_list)
  col_name <- sprintf("dist.p2.ws.g%s.n%d.p%d",
                      format(g, trim = TRUE), N_GRID, p)
  columns[[col_name]] <- xs
  cat(sprintf("    q95=%.3f  cell=%s  total=%s\n\n",
              quantile(xs, 0.95),
              fmt_duration(as.numeric(difftime(Sys.time(), t_cell, units = "secs"))),
              fmt_duration(as.numeric(difftime(Sys.time(), t0,     units = "secs")))))
}

df <- as.data.frame(columns, check.names = FALSE)
if (is.null(OUT_CSV)) {
  OUT_CSV <- file.path(root, "results",
                       sprintf("dist_p2_ws_n%d_t.csv", N_GRID))
}
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("wrote %s (%d rows x %d cols)\n", OUT_CSV, nrow(df), ncol(df)))
