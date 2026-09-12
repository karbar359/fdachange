# Gamma-independent core of the Wilcoxon-type scan. Returns the unweighted
# |H(a, b)| over every window {a..b}, 1 <= a <= b <= n, together with the
# window lengths, the sigma normaliser and the (a, b) index vectors. The
# full window (1, n) has an empty complement, so its increment is 0 and its
# weight is 0; it drops out on its own.
.wilcoxon_core <- function(X, kernel) {
  feature <- switch(kernel,
    h_w = identity, h_0 = abs,
    stop("unknown wilcoxon kernel: ", kernel))
  X <- feature(X)
  n <- length(X)
  H <- outer(X, X, FUN = function(x, y) sign(x - y))
  sigma <- sqrt(mean(colMeans(H)^2))
  if (!is.finite(sigma) || sigma == 0) return(NULL)

  # P[i+1, j+1] = sum over rows 1..i, cols 1..j (plain 1-based summed-area).
  P <- matrix(0, n + 1L, n + 1L)
  P[-1L, -1L] <- t(apply(apply(H, 2L, cumsum), 1L, cumsum))

  b_vec <- rep.int(seq_len(n), seq_len(n))
  a_vec <- sequence(seq_len(n))
  rows  <- P[cbind(b_vec + 1L, n + 1L)] - P[cbind(a_vec, n + 1L)]
  block <- P[cbind(b_vec + 1L, b_vec + 1L)] - P[cbind(b_vec + 1L, a_vec)] -
           P[cbind(a_vec, b_vec + 1L)] + P[cbind(a_vec, a_vec)]
  list(diff_abs = abs(rows - block), l_vec = b_vec - a_vec + 1L,
       a_vec = a_vec, b_vec = b_vec, sigma = sigma, n = n)
}

#' Wilcoxon-type (antisymmetric kernel) scan statistic
#'
#' For each window \eqn{I(k, m) = \{k+1, \dots, m\}} computes
#' \deqn{H(k, m) = \sum_{i \in I, j \notin I} \mathrm{sign}(f(X_i) - f(X_j))}
#' where `f` is the feature map selected by `kernel` (`"h_w"`: identity,
#' sensitive to location shifts; `"h_0"`: absolute value, sensitive to
#' scale/tail changes), and returns the weighted scan
#' \deqn{T = \max_{k < m} |H(k, m)| \, w\big(\tfrac{l}{n}(1 - \tfrac{l}{n})\big)
#'   \, / \, (n^{3/2} \hat\sigma), \qquad l = m - k,}
#' with the symmetric Hoelder weight of [holder_weight()] and
#' \eqn{\hat\sigma^2 = \mathrm{mean}_j (\mathrm{mean}_i \,
#' \mathrm{sign}(f(X_i) - f(X_j)))^2}. All window lengths
#' `1 <= l <= n` are scanned. \eqn{\hat\sigma} is feature-invariant for
#' continuous symmetric distributions, so the same null distribution
#' ([wilcoxon_null_draws()]) applies to both kernel choices.
#'
#' `gamma` may be a vector: the gamma-independent work (the sign matrix,
#' its summed-area table, \eqn{\hat\sigma}) is shared and only the
#' weighting is redone per entry, which is much faster than one call per
#' gamma.
#'
#' @param x numeric sample vector (project functional data first, see
#'   [project_l2()] and friends).
#' @param gamma Hoelder weight exponent(s) in `[0, 0.5]`; may be a vector.
#' @param beta logarithmic weight exponent, see [holder_weight()].
#' @param kernel `"h_w"` (identity feature) or `"h_0"` (absolute value).
#' @return `wilcoxon_statistic()`: numeric vector of scan statistics, one
#'   per entry of `gamma`. `wilcoxon_scan()`: a list with elements `stat`,
#'   `k`, `m` (the maximising window as \eqn{I = \{k+1, \dots, m\}}) and
#'   `interval` (`c(k, m) / n`).
#' @references Bartkus (2026), Chapter 3; Bartkus and Rackauskas (2025).
#' @seealso [epidemic_test()] for the calibrated test,
#'   [wilcoxon_null_draws()] and [kuiper_cdf()] for the null.
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' x[41:70] <- x[41:70] + 1          # epidemic mean shift
#' wilcoxon_statistic(x, gamma = c(0, 0.25))
#' wilcoxon_scan(x, gamma = 0.25)
#' @export
wilcoxon_statistic <- function(x, gamma, beta = 0, kernel = c("h_w", "h_0")) {
  kernel <- match.arg(kernel)
  co <- .wilcoxon_core(x, kernel)
  if (is.null(co)) return(rep(0, length(gamma)))
  norm <- co$n^1.5 * co$sigma
  vapply(gamma, function(g) {
    w <- holder_weights(co$n, g, beta, symmetric = TRUE)
    max(co$diff_abs * w[co$l_vec]) / norm
  }, numeric(1))
}

#' @rdname wilcoxon_statistic
#' @export
wilcoxon_scan <- function(x, gamma, beta = 0, kernel = c("h_w", "h_0")) {
  kernel <- match.arg(kernel)
  stopifnot(length(gamma) == 1L)
  co <- .wilcoxon_core(x, kernel)
  if (is.null(co)) {
    return(list(stat = 0, k = NA_integer_, m = NA_integer_,
                interval = c(NA_real_, NA_real_)))
  }
  w     <- holder_weights(co$n, gamma, beta, symmetric = TRUE)
  vals  <- co$diff_abs * w[co$l_vec]
  i_max <- which.max(vals)
  # Report in the thesis's I(k, m) = {k+1..m} convention: k = a - 1, m = b.
  list(
    stat     = vals[i_max] / (co$n^1.5 * co$sigma),
    k        = co$a_vec[i_max] - 1L,
    m        = co$b_vec[i_max],
    interval = c(co$a_vec[i_max] - 1L, co$b_vec[i_max]) / co$n
  )
}

#' Null draws for the Wilcoxon-type scan
#'
#' Draws of the two-sided Hoelder-weighted Brownian-bridge range
#' \eqn{\sup_{s < t} |W(t) - W(s)| \, w(|t-s|(1-|t-s|))} that calibrate
#' [wilcoxon_statistic()]. At `(gamma = 0, beta = 0)` the law is Kuiper's
#' V and the draws are a deterministic quantile grid of the closed form
#' ([kuiper_cdf()]); other cells come from the shipped Monte-Carlo cache
#' (30,000 draws per cell, gammas 0.05 to 0.495 at `beta = 0`, plus
#' `gamma = 0.5, beta = 1`), extendable via [simulate_wilcoxon_null()] +
#' [save_null_draws()].
#'
#' @inheritParams wilcoxon_statistic
#' @param gamma scalar weight exponent.
#' @return numeric vector of null draws.
#' @examples
#' quantile(wilcoxon_null_draws(gamma = 0), 0.95)   # 1.7473 (Kuiper)
#' @export
wilcoxon_null_draws <- function(gamma, beta = 0) {
  stopifnot(length(gamma) == 1L, length(beta) == 1L)
  if (gamma == 0 && beta == 0) return(.kuiper_grid())
  .null_get("wilcoxon", .null_key("wilcoxon", gamma = gamma, beta = beta),
            sprintf("simulate_wilcoxon_null(n_draws, gamma = %s, beta = %s)",
                    gamma, beta))
}

#' @rdname wilcoxon_null_draws
#' @param alpha significance level.
#' @param splits split-aggregation adjustment, see [critical_value()].
#' @return `wilcoxon_critical_value()`: scalar critical value.
#' @export
wilcoxon_critical_value <- function(gamma, beta = 0, alpha = 0.05, splits = 1) {
  critical_value(wilcoxon_null_draws(gamma, beta), alpha, splits)
}
