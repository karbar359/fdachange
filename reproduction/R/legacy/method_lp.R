#' Lp-type scan test (Chapter 4 of the dissertation).
#'
#' Transforms X to uniform via the scenario-specific null CDF F_0 = `cdf`,
#' then maximises a weighted Lp-norm of the empirical process over all
#' change-windows [k, k+l-1]:
#'   S_p(k,l) = int_0^1 | sum_{i in [k,k+l-1]} ( I(F_0(X_i) <= t) - t ) |^p dt
#'   T_rho   = max_{k+l <= n} S_p(k,l)^(1/p) * w(l/n) / sqrt(n)
#' with non-symmetric weight w(l/n) = rho(l/n).
#'
#' @param cdf vectorised null CDF F_0 of the scenario (S1: pnorm, S2: pchisq,
#'   S3/S5: ppareto, S4: pbeta). Default is standard normal so legacy S1
#'   callers keep working.
#'
#' Callers load shared helpers once (weights.R, null_cache.R).

lp_statistic <- function(X, gamma, p = 2, cdf = function(x) pnorm(x, 0, 1),
                         t_grid = 100, ...) {
  n <- length(X)
  F.X <- cdf(X)
  w <- weight_cache(n, gamma)
  ts <- seq(0, 1, length.out = t_grid + 1)

  # (n x (t_grid+1)) matrix of I(F.X_i <= t_j) - t_j
  ind <- outer(F.X, ts, FUN = function(f, t) (f <= t) - t)

  stat <- 0
  for (k in seq_len(n - 1)) {
    # running cumulative sum along i for this k; start with row k, extend to k+l-1
    cs <- ind[k, ]
    for (l in seq_len(n - k)) {
      if (l > 1) cs <- cs + ind[k + l - 1, ]
      sum_k_l <- sum(abs(cs)^p) / t_grid
      val <- sum_k_l^(1 / p) * w[l]
      if (val > stat) stat <- val
    }
  }
  stat / sqrt(n)
}

#' CORRECTED statistic: restore the F_n centring (2026-09 audit, defect L1).
#'
#' `lp_statistic` faithfully implements the thesis's REDUCED statistic
#' (`ch4_lp_statistics/sec_theory.tex:16`, integrand 1{F(X_i) <= t} - t), but
#' the thesis's reduction from the primary F_n-centred statistic to that form
#' is invalid at first order, and the null cache samples the CENTRED limit
#' (`generators/lp_null.R:89` carries the -(t-s)W(u,1) term). On the
#' parametric-CDF path (scalar scenarios S1-S5, cdf = pnorm etc.) the
#' uncentred statistic against the centred null gives empirical size ~0.4 at
#' the nominal 5% level -- visible in the printed scalar Lp tables' a = 0
#' column (0.40-0.46 at gamma = 0). On the ecdf path (F1-F4, application) the
#' rank transform makes the global sum O(1), supplying the centring
#' accidentally, so this correction is a numerical no-op there.
#'
#' The fix subtracts the in-window share of the global sum, (l/n) * colSums,
#' which matches both the primary printed definition and the cache. The scan
#' stays unrestricted in l: the printed quantile Table 4.1 is cell-exact to
#' the unrestricted cache, so the printed l <= n/2 in the formula contradicts
#' the printed table; we side with the table and record the erratum. Length-1
#' windows starting at n are included for completeness.
lp_statistic_corrected <- function(X, gamma, p = 2,
                                   cdf = function(x) pnorm(x, 0, 1),
                                   t_grid = 100, ...) {
  n <- length(X)
  F.X <- cdf(X)
  w <- weight_cache(n, gamma)
  ts <- seq(0, 1, length.out = t_grid + 1)
  ind <- outer(F.X, ts, FUN = function(f, t) (f <= t) - t)
  tot <- colSums(ind)                       # F_n centring vector

  stat <- 0
  for (k in seq_len(n)) {
    cs <- numeric(t_grid + 1)
    for (l in seq_len(n - k + 1)) {
      cs <- cs + ind[k + l - 1, ]
      cc <- cs - (l / n) * tot
      sum_k_l <- sum(abs(cc)^p) / t_grid
      val <- sum_k_l^(1 / p) * w[l]
      if (val > stat) stat <- val
    }
  }
  stat / sqrt(n)
}

#' Corrected scan with argmax: same statistic, additionally returns the
#' maximising window as (k, m) with I = {k+1..m} and the fractional interval.
#' Needed by the application tables -- the legacy pipeline never exposed the
#' Lp window; the printed Table 7.1 endpoints came from the archive script.
lp_scan <- function(X, gamma, p = 2, cdf = function(x) pnorm(x, 0, 1),
                    t_grid = 100, ...) {
  n <- length(X)
  F.X <- cdf(X)
  w <- weight_cache(n, gamma)
  ts <- seq(0, 1, length.out = t_grid + 1)
  ind <- outer(F.X, ts, FUN = function(f, t) (f <= t) - t)
  tot <- colSums(ind)

  stat <- 0; bk <- NA_integer_; bm <- NA_integer_
  for (k in seq_len(n)) {
    cs <- numeric(t_grid + 1)
    for (l in seq_len(n - k + 1)) {
      cs <- cs + ind[k + l - 1, ]
      cc <- cs - (l / n) * tot
      val <- (sum(abs(cc)^p) / t_grid)^(1 / p) * w[l]
      if (val > stat) { stat <- val; bk <- k - 1L; bm <- k + l - 1L }
    }
  }
  list(stat = stat / sqrt(n), k = bk, m = bm, interval = c(bk, bm) / n)
}

#' Null distribution: precomputed Wiener-sheet draws from `.raw/dist_p2_ws_n200_t.csv`.
#' Indexed by (gamma, p); n = 200 fixed in the cache.
lp_null_distribution <- function(n, gamma, p = 2, ...) {
  load_lp_null(gamma, p)
}
