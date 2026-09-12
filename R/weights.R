#' Hoelder-type scan weight
#'
#' The weight applied to a scan window of relative length `t` in `(0, 1)`:
#' \deqn{w(t) = t^{-\gamma} \, \log(c/t)^{-\beta}, \qquad c = 0.3,}
#' with `w(0) = 0` by convention. The scan statistics **multiply** the
#' window increment by this weight, so larger `gamma` up-weights short
#' windows. The logarithmic factor is only active at `beta != 0`, which
#' the dissertation pairs exclusively with `gamma = 0.5`.
#'
#' @param t numeric vector of relative window lengths in `[0, 1)`.
#' @param gamma power-weight exponent, `0 <= gamma <= 0.5`.
#' @param beta logarithmic-weight exponent (default 0).
#' @return numeric vector of weights, same length as `t`.
#' @seealso [holder_weights()] for the per-window-length cache used by the
#'   scan statistics.
#' @examples
#' holder_weight(c(0.1, 0.25, 0.5), gamma = 0.25)
#' @export
holder_weight <- function(t, gamma, beta = 0) {
  denom <- t^gamma * log(0.3 / t)^beta
  ifelse(t != 0, 1 / denom, 0)
}

#' Weight cache over integer window lengths
#'
#' Precomputes [holder_weight()] for every window length `l = 1, ..., n - 1`
#' of a sample of size `n`; entry `n` is 0 (a full-sample window has an
#' empty complement and drops out of every scan).
#'
#' @param n sample size.
#' @param gamma,beta weight exponents, see [holder_weight()].
#' @param symmetric if `FALSE` (default), `w[l] = w(l/n)`: the direct
#'   weight used by the RKHS and Lp scans. If `TRUE`,
#'   `w[l] = w((l/n)(1 - l/n))`: the symmetrised weight used by the
#'   Wilcoxon-type scan.
#' @return numeric vector of length `n`; `w[l]` is the weight for a
#'   window holding `l` observations.
#' @examples
#' w <- holder_weights(100, gamma = 0.25, symmetric = TRUE)
#' plot(w, type = "l", xlab = "window length l", ylab = "weight")
#' @export
holder_weights <- function(n, gamma, beta = 0, symmetric = FALSE) {
  w <- numeric(n)
  for (l in seq_len(n - 1)) {
    t <- l / n
    if (symmetric) t <- t * (1 - t)
    w[l] <- holder_weight(t, gamma, beta)
  }
  w
}

#' Critical value from Monte-Carlo null draws
#'
#' Upper quantile of a null sample, with the split-aggregation adjustment
#' of [scan_over_splits()]: for a statistic that is the maximum over
#' `splits` (asymptotically independent) chunk statistics, the level-alpha
#' critical value is the `(1 - alpha)^(1/splits)` quantile of the
#' single-chunk null.
#'
#' The quantile is computed with [stats::quantile()]'s default type 7;
#' published thresholds depend on that choice, so it is deliberately not
#' configurable.
#'
#' @param draws numeric vector of null-distribution draws.
#' @param alpha significance level (default 0.05).
#' @param splits number of sample splits the statistic aggregates over
#'   (default 1, no adjustment).
#' @return scalar critical value.
#' @examples
#' critical_value(wilcoxon_null_draws(gamma = 0), alpha = 0.05)
#' @export
critical_value <- function(draws, alpha = 0.05, splits = 1) {
  as.numeric(quantile(draws, probs = (1 - alpha)^(1 / splits)))
}
