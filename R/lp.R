#' Lp-type scan statistic
#'
#' Transforms the sample to uniform via the null CDF `cdf`, then maximises
#' a weighted Lp norm of the centred empirical process over all change
#' windows \eqn{I = \{k+1, \dots, k+l\}}:
#' \deqn{S_p(k, l) = \int_0^1 \Big| \sum_{i \in I} \big(\mathbf{1}\{F_0(X_i)
#'   \le t\} - t\big) - \tfrac{l}{n} \sum_{i=1}^n \big(\mathbf{1}\{F_0(X_i)
#'   \le t\} - t\big) \Big|^p dt,}
#' \deqn{T = \max_{k, l} S_p(k, l)^{1/p} \, w(l/n) / \sqrt{n},}
#' with the non-symmetric Hoelder weight of [holder_weight()]. The
#' \eqn{F_n}-centring term (the second sum) makes the statistic match the
#' Wiener-sheet null of [lp_null_draws()]; on the ecdf path it is a
#' numerical no-op, on a parametric-CDF path it is essential.
#'
#' With `cdf = NULL` (default) the empirical CDF of the sample itself is
#' used: the safe choice for real data and for projected functional
#' scenarios. Pass the true null CDF when it is known (e.g.
#' `scenario_cdf("S1")` for a standard normal null).
#'
#' @param x numeric sample vector.
#' @param gamma Hoelder weight exponent (scalar).
#' @param p Lp exponent (the shipped null covers p = 2, 4, 10).
#' @param cdf vectorised null CDF `F_0`, or `NULL` for the empirical CDF.
#' @param t_grid number of integration subintervals on `[0, 1]`.
#' @return `lp_statistic()`: scalar statistic. `lp_scan()`: a list with
#'   `stat`, `k`, `m` (maximising window \eqn{I = \{k+1, \dots, m\}}) and
#'   `interval` (`c(k, m) / n`).
#' @references Bartkus (2026), Chapter 4; Bartkus and Rackauskas (2026).
#' @seealso [epidemic_test()], [lp_null_draws()], [scenario_cdf()]
#' @examples
#' set.seed(1)
#' x <- rnorm(120)
#' x[49:84] <- x[49:84] + 1
#' lp_statistic(x, gamma = 0.25, cdf = pnorm)
#' lp_scan(x, gamma = 0.25)         # ecdf path
#' @export
lp_statistic <- function(x, gamma, p = 2, cdf = NULL, t_grid = 100) {
  n <- length(x)
  F.X <- if (is.null(cdf)) ecdf(x)(x) else cdf(x)
  w <- holder_weights(n, gamma)
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

#' @rdname lp_statistic
#' @export
lp_scan <- function(x, gamma, p = 2, cdf = NULL, t_grid = 100) {
  n <- length(x)
  F.X <- if (is.null(cdf)) ecdf(x)(x) else cdf(x)
  w <- holder_weights(n, gamma)
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

#' Null draws for the Lp-type scan
#'
#' Monte-Carlo draws of the weighted Lp functional of the centred Wiener
#' sheet that calibrates [lp_statistic()]. The shipped cache holds 300
#' draws per cell on a 200-point grid for `gamma` in \{0, 0.25, 0.45\}
#' and `p` in \{2, 4, 10\}; extend via [simulate_lp_null()] +
#' [save_null_draws()].
#'
#' @param gamma Hoelder weight exponent (scalar).
#' @param p Lp exponent.
#' @return numeric vector of null draws.
#' @examples
#' quantile(lp_null_draws(gamma = 0.25, p = 2), 0.95)
#' @export
lp_null_draws <- function(gamma, p = 2) {
  stopifnot(length(gamma) == 1L, length(p) == 1L)
  .null_get("lp", .null_key("lp", gamma = gamma, p = p),
            sprintf("simulate_lp_null(n_draws, gamma = %s, p = %s)", gamma, p))
}

#' @rdname lp_null_draws
#' @param alpha significance level.
#' @param splits split-aggregation adjustment, see [critical_value()].
#' @return `lp_critical_value()`: scalar critical value.
#' @export
lp_critical_value <- function(gamma, p = 2, alpha = 0.05, splits = 1) {
  critical_value(lp_null_draws(gamma, p), alpha, splits)
}
