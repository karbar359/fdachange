#' Rho weighting on (0,1): rho(t) = 1 / ( t^gamma * log(c/t)^beta ).
rho_gamma_beta <- function(t, gamma, beta = 0) {
  c.w <- .25 + (.5 - 0.45)
  denom <- t^gamma * log(c.w / t)^beta
  ifelse(t != 0, 1 / denom, 0)
}

#' Cache weights for a scan statistic over windows of size l = 1,...,n-1.
#'
#' @param symmetric if FALSE (default), w[l] = rho(l/n) -- the direct
#'   Hoelder-type weight used by the RKHS and L_p tests (Ch4, Ch5).
#'   If TRUE, w[l] = rho(l(n-l)/n^2) -- the symmetrised weight used by
#'   the antisymmetric-kernel Wilcoxon test (Ch3).
weight_cache <- function(n, gamma, beta = 0, symmetric = FALSE) {
  w <- numeric(n)
  for (l in seq_len(n - 1)) {
    t <- l / n
    if (symmetric) t <- t * (1 - t)
    w[l] <- rho_gamma_beta(t, gamma, beta)
  }
  w
}

#' Upper alpha-quantile of a Monte-Carlo null distribution.
critical_value <- function(null_dist, alpha = 0.05, m.d = 1) {
  as.numeric(quantile(null_dist, probs = (1 - alpha)^(1 / m.d)))
}
