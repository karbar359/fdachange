#' Kuiper law: the exact gamma = 0 Wilcoxon-scan null
#'
#' At `gamma = 0` the Hoelder weight is identically 1 and the Wilcoxon
#' scan statistic converges to the range of a standard Brownian bridge,
#' whose distribution is Kuiper's V:
#' \deqn{P(V \le v) = 1 - \sum_{k \ge 1} 2\,(4 k^2 v^2 - 1)\, e^{-2 k^2 v^2}.}
#' `kuiper_cdf()` evaluates that series; `kuiper_quantile()` inverts it by
#' root finding. The package always serves the `(gamma = 0, beta = 0)`
#' Wilcoxon null from this closed form rather than from a Monte-Carlo
#' cache.
#'
#' @param v numeric vector of quantiles.
#' @param p numeric vector of probabilities.
#' @return `kuiper_cdf()`: vector of cumulative probabilities.
#'   `kuiper_quantile()`: vector of quantiles.
#' @examples
#' kuiper_quantile(0.95)   # 1.7473
#' kuiper_cdf(1.7473)
#' @export
kuiper_cdf <- function(v) {
  vapply(v, function(vi) {
    if (vi <= 0.3) return(0)
    k <- 1:1000
    min(1, max(0, 1 - sum(2 * (4 * k^2 * vi^2 - 1) * exp(-2 * k^2 * vi^2))))
  }, numeric(1))
}

#' @rdname kuiper_cdf
#' @export
kuiper_quantile <- function(p) {
  vapply(p, function(pp) {
    pp <- min(max(pp, 1e-9), 1 - 1e-12)
    uniroot(function(v) kuiper_cdf(v) - pp, c(0.3, 6),
            tol = .Machine$double.eps^0.5)$root
  }, numeric(1))
}

# Deterministic quantile grid standing in for a Monte-Carlo sample of the
# Kuiper law: qK((i - 0.5)/N), i = 1..N. Empirical quantiles of this vector
# reproduce the exact quantiles to O(1/N), so every consumer of
# wilcoxon_null_draws() (which all take empirical quantiles) works
# unchanged. Built by inverting a dense tabulated CDF (vectorised, ~0.1 s)
# rather than 200k uniroot calls; interpolation error is ~1e-5, far below
# the 1e-3 the 3-decimal critical values need. kuiper_quantile() stays
# uniroot-exact for spot checks.
.kuiper_env <- new.env(parent = emptyenv())
.kuiper_grid <- function(N = 200000L) {
  if (is.null(.kuiper_env$g)) {
    v   <- seq(0.3, 6, length.out = 20001L)
    vk  <- outer(v^2, (1:100)^2)          # v^2 * k^2 grid
    cdf <- 1 - rowSums(2 * (4 * vk - 1) * exp(-2 * vk))
    cdf <- cummax(pmin(1, pmax(0, cdf)))  # clamp + enforce monotonicity
    keep <- !duplicated(cdf)
    .kuiper_env$g <- approx(cdf[keep], v[keep],
                            xout = (seq_len(N) - 0.5) / N, rule = 2)$y
  }
  .kuiper_env$g
}
