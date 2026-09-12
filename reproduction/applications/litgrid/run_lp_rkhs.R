#!/usr/bin/env Rscript
# L_p and RKHS scans on the Litgrid balancing-price data, both with and
# without the outlier week. Scalar input is the L_2-norm projection of
# each weekly curve (matches sec_framework.tex's "the observed sample
# X_1,..,X_n"); the F_n / chi^2 transforms of T1, T2 are then applied
# before passing to the test statistic.
#
# Output: TSV-like rows that can be pasted into tab:litgrid_asym
# (with-outlier) and used to re-validate tab:litgrid-combined
# (no-outlier).
#
# Usage:
#   Rscript applications/litgrid/run_lp_rkhs.R

ALPHA      <- 0.05
ROW_DROP   <- 204L            # outlier week (excluded in "no outlier" pass)
ROW_KEEP_N <- 344L
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

# --- Local scan helpers (track argmax, unlike methods/{lp,rkhs}/run.R which
#     return only the test statistic) ----------------------------------------

lp_scan <- function(X, gamma, p = 2L, t_grid = 100L) {
  n   <- length(X)
  F.X <- X            # caller passes already-transformed (F_n) ranks in [0,1]
  w   <- weight_cache(n, gamma)
  ts  <- seq(0, 1, length.out = t_grid + 1L)
  ind <- outer(F.X, ts, FUN = function(f, t) (f <= t) - t)
  best_stat <- 0; best_k <- NA_integer_; best_l <- NA_integer_
  for (k in seq_len(n - 1L)) {
    cs <- ind[k, ]
    for (l in seq_len(n - k)) {
      if (l > 1L) cs <- cs + ind[k + l - 1L, ]
      sum_k_l <- sum(abs(cs)^p) / t_grid
      val     <- sum_k_l^(1 / p) * w[l]
      if (val > best_stat) { best_stat <- val; best_k <- k; best_l <- l }
    }
  }
  list(stat = best_stat / sqrt(n),
       k = best_k, m = best_k + best_l - 1L,
       interval = c((best_k - 1L) / n, (best_k + best_l - 1L) / n))
}

rkhs_scan <- function(X, gamma, kernel_name = "gaussian") {
  n  <- length(X)
  K  <- gram(X, kernel_fn(kernel_name))
  K.sum <- sum(K)
  sat   <- NumMatrix(K)$sat
  w     <- weight_cache(n, gamma)
  n1    <- n - 1L
  m_vec <- rep(seq_len(n1), times = seq_len(n1))
  k_vec <- sequence(seq_len(n1))
  l_vec <- m_vec - k_vec + 1L
  r2    <- pmin(m_vec + 2L, nrow(sat))
  c2    <- pmin(m_vec + 2L, ncol(sat))
  k1    <- k_vec + 1L
  c_end <- rep(ncol(sat), length(k_vec))
  two   <- rep(2L,         length(k_vec))
  cross <- sat[cbind(r2, c_end)] - sat[cbind(r2, two)] -
           sat[cbind(k1, c_end)] + sat[cbind(k1, two)]
  B_km  <- sat[cbind(r2, c2)]   - sat[cbind(r2, k1)]   -
           sat[cbind(k1, c2)]   + sat[cbind(k1, k1)]
  l_n   <- l_vec / n
  s_km  <- B_km + l_n^2 * K.sum - 2 * l_n * cross
  s_km[s_km < 0] <- 0
  vals  <- sqrt(s_km) * w[l_vec]
  idx   <- which.max(vals)
  list(stat = vals[idx] / sqrt(n),
       k = k_vec[idx], m = m_vec[idx],
       interval = c((k_vec[idx] - 1L) / n, m_vec[idx] / n))
}

# --- Data -----------------------------------------------------------------

litgrid_dir <- Sys.getenv("LITGRID_RAW_DIR", unset = "")
if (!nzchar(litgrid_dir)) {
  litgrid_dir <- normalizePath(file.path(root, "..", "data-raw", "litgrid-raw"),
                               mustWork = FALSE)
}
if (!dir.exists(litgrid_dir)) {
  stop("Litgrid raw CSVs not found at ", litgrid_dir,
       "; run Rscript applications/litgrid/download_litgrid.R or set LITGRID_RAW_DIR.")
}
cat(sprintf("Loading Litgrid balancing prices from %s\n", litgrid_dir))
f_matrix_full <- load_litgrid_balancing_prices(litgrid_dir)
cat(sprintf("Raw weekly matrix: %d weeks x %d hours\n",
            nrow(f_matrix_full), ncol(f_matrix_full)))


run_pass <- function(label, f_matrix) {
  X <- t(f_matrix)
  n <- ncol(X)
  set.seed(SEED)
  X_proj <- projection_fn("l2")(X)
  cat(sprintf("\n=== %s pass: n=%d, l2-norm projection (sd=%.3g) ===\n",
              label, n, sd(X_proj)))

  ranks  <- rank(X_proj, ties.method = "average") / n      # T1: F_n transform
  chi_in <- qchisq(pmin(ranks, 1 - 1e-9), df = 1)          # T2: chi^2_1 quantile

  cat("method\tgamma\tT\tk_idx\tm_idx\tlow\thigh\tcv95\treject\n")
  for (gamma in c(0, 0.25, 0.45)) {
    sl <- lp_scan(ranks, gamma, p = 2L)
    nl <- lp_null_distribution(n, gamma, p = 2L)
    cvl <- critical_value(nl, ALPHA, m.d = 1L)
    cat(sprintf("lp\t%.2f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%s\n",
                gamma, sl$stat, sl$k, sl$m, sl$interval[1], sl$interval[2],
                cvl, sl$stat > cvl))
  }
  for (gamma in c(0, 0.25, 0.45)) {
    sr <- rkhs_scan(chi_in, gamma, kernel_name = "gaussian")
    nr <- rkhs_null_distribution(n, gamma, kernel_name = "gaussian")
    cvr <- critical_value(nr, ALPHA, m.d = 1L)
    cat(sprintf("rkhs\t%.2f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%s\n",
                gamma, sr$stat, sr$k, sr$m, sr$interval[1], sr$interval[2],
                cvr, sr$stat > cvr))
  }
}

# No-outlier pass: drop row 204, keep first 344 of remainder.
kept <- setdiff(seq_len(nrow(f_matrix_full)), ROW_DROP)
m_no <- f_matrix_full[kept, , drop = FALSE][seq_len(ROW_KEEP_N), , drop = FALSE]
run_pass("NO OUTLIER (drop row 204)", m_no)

# With-outlier pass: keep first 344 rows as-is.
m_with <- f_matrix_full[seq_len(ROW_KEEP_N), , drop = FALSE]
run_pass("WITH OUTLIER (full dataset)", m_with)
