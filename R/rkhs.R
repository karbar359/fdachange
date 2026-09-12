# Core of the RKHS scan: weighted embedding distances over all windows
# I = {a..b} of length l = b - a + 1 <= n/2, all starts a >= 1.
.rkhs_core <- function(X, gamma, kernel, beta = 0) {
  n <- length(X)
  if (n < 2L) return(NULL)
  K <- gram(X, .resolve_kernel(kernel))
  P <- matrix(0, n + 1L, n + 1L)
  P[-1L, -1L] <- t(apply(apply(K, 2L, cumsum), 1L, cumsum))
  K.sum <- P[n + 1L, n + 1L]
  w <- holder_weights(n, gamma, beta)

  l_max <- n %/% 2L
  starts_per_l <- n - seq_len(l_max) + 1L         # l has n - l + 1 start points
  a_vec <- sequence(starts_per_l)
  l_rep <- rep.int(seq_len(l_max), starts_per_l)
  b_vec <- a_vec + l_rep - 1L
  rows  <- P[cbind(b_vec + 1L, n + 1L)] - P[cbind(a_vec, n + 1L)]
  block <- P[cbind(b_vec + 1L, b_vec + 1L)] - P[cbind(b_vec + 1L, a_vec)] -
           P[cbind(a_vec, b_vec + 1L)] + P[cbind(a_vec, a_vec)]
  l_n   <- l_rep / n
  s_km  <- block + l_n^2 * K.sum - 2 * l_n * rows
  s_km[s_km < 0] <- 0                     # clamp floating-point negatives
  list(vals = sqrt(s_km) * w[l_rep], a_vec = a_vec, b_vec = b_vec, n = n)
}

#' RKHS (kernel mean embedding) scan statistic
#'
#' For each window \eqn{I = \{a, \dots, b\}} of length
#' \eqn{l = b - a + 1 \le n/2} computes the squared embedding distance
#' \deqn{\Delta^2(I) = \sum_{i, j \in I} K_{ij}
#'   - 2 \tfrac{l}{n} \sum_{i \in I}\sum_{j=1}^n K_{ij}
#'   + \big(\tfrac{l}{n}\big)^2 \sum_{i, j = 1}^n K_{ij}}
#' and returns the weighted scan
#' \deqn{T = \max_{I} \sqrt{\Delta^2(I)} \; w(l/n) / \sqrt{n}}
#' with the non-symmetric Hoelder weight of [holder_weight()].
#'
#' @param x numeric sample vector (project functional data first).
#' @param gamma Hoelder weight exponent (scalar).
#' @param kernel a kernel name (`"gaussian"`, `"min"`, `"2max"`) or a
#'   kernel function of two vectorised arguments.
#' @param beta logarithmic weight exponent, see [holder_weight()].
#' @return `rkhs_statistic()`: scalar statistic. `rkhs_scan()`: a list
#'   with `stat`, `k`, `m` (maximising window \eqn{I = \{k+1, \dots, m\}})
#'   and `interval` (`c(k, m) / n`).
#' @references Bartkus (2026), Chapter 5; Bartkus and Rackauskas (2026).
#' @seealso [epidemic_test()], [rkhs_null_draws()], [kernel_gaussian()]
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' x[41:70] <- x[41:70] + 1
#' rkhs_statistic(x, gamma = 0.25, kernel = "gaussian")
#' rkhs_scan(x, gamma = 0.25, kernel = "gaussian")
#' @export
rkhs_statistic <- function(x, gamma, kernel = "gaussian", beta = 0) {
  co <- .rkhs_core(x, gamma, kernel, beta)
  if (is.null(co)) return(0)
  max(co$vals) / sqrt(co$n)
}

#' @rdname rkhs_statistic
#' @export
rkhs_scan <- function(x, gamma, kernel = "gaussian", beta = 0) {
  co <- .rkhs_core(x, gamma, kernel, beta)
  if (is.null(co)) {
    return(list(stat = 0, k = NA_integer_, m = NA_integer_,
                interval = c(NA_real_, NA_real_)))
  }
  i_max <- which.max(co$vals)
  list(stat     = co$vals[i_max] / sqrt(co$n),
       k        = co$a_vec[i_max] - 1L,
       m        = co$b_vec[i_max],
       interval = c(co$a_vec[i_max] - 1L, co$b_vec[i_max]) / co$n)
}

#' Null draws for the RKHS scan
#'
#' Monte-Carlo draws of [rkhs_statistic()] under an i.i.d. null at
#' `n = 2000`, keyed by `(gamma, kernel, baseline)`. The `baseline` tag
#' names the H0 sampling distribution of the cached cell: `""` is
#' standard normal (shipped for the Gaussian kernel), `"runif"` is
#' U(0, 1) (shipped for min and 2max), `"rnorm.sq"` is a squared standard
#' normal (min), and `"rbeta.1.9"` is Beta(1, 9) (2max). Shipped gammas:
#' 0, 0.25, 0.45. Other baselines or gammas: [simulate_rkhs_null()] +
#' [save_null_draws()].
#'
#' @param gamma Hoelder weight exponent (scalar).
#' @param kernel kernel name: `"gaussian"`, `"min"` or `"2max"`.
#' @param baseline H0-distribution tag, see Details.
#' @return numeric vector of null draws.
#' @examples
#' quantile(rkhs_null_draws(gamma = 0.25, kernel = "gaussian"), 0.95)
#' @export
rkhs_null_draws <- function(gamma, kernel = "gaussian", baseline = "") {
  stopifnot(length(gamma) == 1L)
  key <- .null_key("rkhs", gamma = gamma, kernel = kernel, baseline = baseline)
  .null_get("rkhs", key,
            sprintf(paste0("simulate_rkhs_null(n_draws, gamma = %s, kernel = ",
                           "\"%s\", rdist = <H0 sampler>)"),
                    gamma, .kernel_name(kernel)))
}

#' @rdname rkhs_null_draws
#' @param alpha significance level.
#' @param splits split-aggregation adjustment, see [critical_value()].
#' @return `rkhs_critical_value()`: scalar critical value.
#' @export
rkhs_critical_value <- function(gamma, kernel = "gaussian", baseline = "",
                                alpha = 0.05, splits = 1) {
  critical_value(rkhs_null_draws(gamma, kernel, baseline), alpha, splits)
}
