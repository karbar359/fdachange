#!/usr/bin/env Rscript
# L_p and RKHS scans on the Litgrid balancing-price data, mirroring the
# Wilcoxon coverage: every test runs over the three projections {ito,
# l2, c_max_var} for both the no-outlier and with-outlier passes.
#
# Projection-kernel pairings (matches F4 simulation calibration):
#   ito       + Gaussian   T1 (rank/n) -> T2 (qnorm)
#   l2        + min        T1 (rank/n) only
#   c_max_var + Gaussian   T1 (rank/n) -> T2 (qnorm)
#
# L_p uses T1 only for all three projections (its asymptotic null is
# parameter-free under U(0,1)).
#
# The ito projection is stochastic in the random-feature draw, so we
# repeat ITO_REPS = 10 times and report mean +/- s.d. of statistic,
# interval endpoints, and rejection rate -- matching the Wilcoxon "It^o
# integral" rows already in tab:litgrid-combined / tab:litgrid_asym.

ALPHA      <- 0.05
ROW_DROP   <- 204L
ROW_KEEP_N <- 344L
SEED       <- 3L
ITO_REPS   <- 10L

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

# --- Local scan helpers (track argmax) -------------------------------------

lp_scan <- function(X, gamma, p = 2L, t_grid = 100L) {
  n   <- length(X)
  F.X <- X
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

# --- Transforms ------------------------------------------------------------

T1_uniform <- function(x) rank(x, ties.method = "average") / length(x)
T2_normal  <- function(x) {
  u <- pmin(pmax(T1_uniform(x), 1e-9), 1 - 1e-9)
  qnorm(u)
}
T2_chi2    <- function(x) {
  u <- pmin(pmax(T1_uniform(x), 1e-9), 1 - 1e-9)
  qchisq(u, df = 1)
}

# --- Single-projection runner ----------------------------------------------

run_one <- function(X_proj, gamma, kernel_for_rkhs, transform_for_rkhs,
                    rkhs_baseline = "") {
  n <- length(X_proj)
  ranks_uniform <- T1_uniform(X_proj)
  rkhs_input    <- transform_for_rkhs(X_proj)

  sl  <- lp_scan(ranks_uniform, gamma, p = 2L)
  nl  <- lp_null_distribution(n, gamma, p = 2L)
  cvl <- critical_value(nl, ALPHA, m.d = 1L)

  sr  <- rkhs_scan(rkhs_input, gamma, kernel_name = kernel_for_rkhs)
  nr  <- rkhs_null_distribution(n, gamma, kernel_name = kernel_for_rkhs,
                                baseline = rkhs_baseline)
  cvr <- critical_value(nr, ALPHA, m.d = 1L)

  list(
    lp = list(stat = sl$stat, low = sl$interval[1], high = sl$interval[2],
              cv = cvl, reject = sl$stat > cvl),
    rkhs = list(stat = sr$stat, low = sr$interval[1], high = sr$interval[2],
                cv = cvr, reject = sr$stat > cvr)
  )
}

# --- Full pass over the three projections ---------------------------------

run_pass <- function(label, f_matrix) {
  X <- t(f_matrix)
  n <- ncol(X)
  cat(sprintf("\n############################\n# %s pass: n=%d\n############################\n",
              label, n))

  # l2 (deterministic): min kernel, T1 only
  set.seed(SEED)
  X_l2 <- projection_fn("l2")(X)
  cat(sprintf("\n=== l2 projection (sd=%.3g) ===\n", sd(X_l2)))
  cat("method\tprojection\tgamma\tT\tlow\thigh\tcv95\treject\n")
  for (gamma in c(0, 0.25, 0.45)) {
    res <- run_one(X_l2, gamma, kernel_for_rkhs = "min",
                   transform_for_rkhs = T2_chi2,
                   rkhs_baseline = "rnorm.sq")
    cat(sprintf("lp\tl2\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\n",
                gamma, res$lp$stat, res$lp$low, res$lp$high, res$lp$cv, res$lp$reject))
    cat(sprintf("rkhs\tl2\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\n",
                gamma, res$rkhs$stat, res$rkhs$low, res$rkhs$high, res$rkhs$cv, res$rkhs$reject))
  }

  # c_max_var (single random projection point, deterministic seed)
  set.seed(SEED)
  X_cmv <- projection_fn("c_max_var")(X)
  cat(sprintf("\n=== c_max_var projection (sd=%.3g) ===\n", sd(X_cmv)))
  cat("method\tprojection\tgamma\tT\tlow\thigh\tcv95\treject\n")
  for (gamma in c(0, 0.25, 0.45)) {
    res <- run_one(X_cmv, gamma, kernel_for_rkhs = "gaussian",
                   transform_for_rkhs = T2_normal)
    cat(sprintf("lp\tc_max_var\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\n",
                gamma, res$lp$stat, res$lp$low, res$lp$high, res$lp$cv, res$lp$reject))
    cat(sprintf("rkhs\tc_max_var\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\n",
                gamma, res$rkhs$stat, res$rkhs$low, res$rkhs$high, res$rkhs$cv, res$rkhs$reject))
  }

  # ito (10 random-feature draws -> mean +/- sd, rejection rate)
  cat(sprintf("\n=== ito projection (10 draws) ===\n"))
  cat("method\tprojection\tgamma\tT_mean\tT_sd\tlow_mean\tlow_sd\thigh_mean\thigh_sd\treject_pct\n")
  set.seed(SEED)
  draws <- vector("list", ITO_REPS)
  for (r in seq_len(ITO_REPS)) {
    draws[[r]] <- projection_fn("ito")(X)
  }
  for (gamma in c(0, 0.25, 0.45)) {
    lp_T <- numeric(ITO_REPS); lp_lo <- lp_T; lp_hi <- lp_T; lp_rej <- logical(ITO_REPS)
    rk_T <- numeric(ITO_REPS); rk_lo <- rk_T; rk_hi <- rk_T; rk_rej <- logical(ITO_REPS)
    for (r in seq_len(ITO_REPS)) {
      res <- run_one(draws[[r]], gamma, kernel_for_rkhs = "gaussian",
                     transform_for_rkhs = T2_normal)
      lp_T[r]  <- res$lp$stat;  lp_lo[r] <- res$lp$low;  lp_hi[r] <- res$lp$high; lp_rej[r] <- res$lp$reject
      rk_T[r]  <- res$rkhs$stat; rk_lo[r] <- res$rkhs$low; rk_hi[r] <- res$rkhs$high; rk_rej[r] <- res$rkhs$reject
    }
    cat(sprintf("lp\tito\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d%%\n",
                gamma, mean(lp_T), sd(lp_T), mean(lp_lo), sd(lp_lo),
                mean(lp_hi), sd(lp_hi), sum(lp_rej) * 10L))
    cat(sprintf("rkhs\tito\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d%%\n",
                gamma, mean(rk_T), sd(rk_T), mean(rk_lo), sd(rk_lo),
                mean(rk_hi), sd(rk_hi), sum(rk_rej) * 10L))
  }
}

# --- Data ------------------------------------------------------------------

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

kept  <- setdiff(seq_len(nrow(f_matrix_full)), ROW_DROP)
m_no  <- f_matrix_full[kept, , drop = FALSE][seq_len(ROW_KEEP_N), , drop = FALSE]
m_with <- f_matrix_full[seq_len(ROW_KEEP_N), , drop = FALSE]

run_pass("NO OUTLIER (drop row 204)", m_no)
run_pass("WITH OUTLIER (full dataset)", m_with)
