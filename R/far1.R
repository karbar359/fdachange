#' Birbilas's functional epidemic-change statistic
#'
#' The CUSUM-type statistic of Birbilas (2025) for functional
#' observations:
#' \deqn{E_{n,\beta,\gamma} = \max_{0 \le k < m \le n} \; \max_{t \in [0,1]}
#'   \frac{\big| (I - \beta(t)) \sum_{i=k+1}^m (X_i(t) - \bar X_n(t)) \big|}
#'        {(m - k)^\gamma},}
#' returned on the test scale \eqn{n^{-1/2+\gamma} E_{n,\beta,\gamma}}.
#' `beta_t` is the pointwise multiplication operator \eqn{\beta(t)};
#' Birbilas's canonical choice is `beta_t = far1_betas$beta1`
#' (\eqn{0.7\sqrt{t}}) with `gamma = 0.4`. Reject for large values against
#' the asymptotic null of [far1_null_draws()], or against a
#' finite-sample null passed to [epidemic_test()] via `null=`.
#'
#' Unlike the other three statistics, this one consumes the functional
#' sample directly (no projection): a numeric matrix with one observation
#' per **column**, rows being the evaluation grid on `[0, 1]`.
#'
#' @param X numeric matrix, `(grid points) x n`, one curve per column.
#' @param gamma window-length exponent; Birbilas's default is 0.4.
#' @param beta_t function of `t` on `[0, 1]`: the operator
#'   \eqn{\beta(t)}. See [far1_betas].
#' @return `far1_statistic()`: scalar statistic (test scale).
#'   `far1_scan()`: a list with `stat`, `k`, `m` (maximising window
#'   \eqn{I = \{k+1, \dots, m\}} in observation indices) and `interval`
#'   (`c(k, m) / n`).
#' @references Birbilas, A. (2025). *Functional Data Analysis in
#'   Telecommunications*. PhD dissertation, Vilnius University.
#'   \doi{10.15388/vu.thesis.841}. Bartkus (2026), Chapter 6.
#' @seealso [epidemic_test()], [far1_null_draws()], [generate_f4()]
#' @examples
#' set.seed(1)
#' X <- generate_f4(60, a = 2, ec_interval = c(0.25, 0.45), l = 50)
#' far1_statistic(X, gamma = 0.4)
#' far1_scan(X, gamma = 0.4)
#' @export
far1_statistic <- function(X, gamma = 0.4,
                               beta_t = far1_betas$beta1) {
  if (!is.matrix(X)) stop("far1_statistic expects a matrix (grid x n)")
  L <- nrow(X); n <- ncol(X)
  if (n < 1L) return(0)
  s <- seq(0, 1, length.out = L)
  W <- (1 - beta_t(s)) * (X - rowMeans(X))   # (I - beta) on centred curves

  # Column prefix sums per row: CS[, j + 1] = sum_{i=1..j} W[, i].
  CS <- cbind(0, t(apply(W, 1L, cumsum)))    # L x (n + 1)

  E_max <- 0
  for (d in 0:(n - 1L)) {                     # window of d + 1 columns
    m_idx <- seq.int(1L, n - d)
    k_idx <- m_idx + d
    block <- abs(CS[, k_idx + 1L, drop = FALSE] - CS[, m_idx, drop = FALSE])
    val_d <- max(block) / (d + 1L)^gamma
    if (val_d > E_max) E_max <- val_d
  }
  E_max * n^(-0.5 + gamma)
}

#' @rdname far1_statistic
#' @export
far1_scan <- function(X, gamma = 0.4, beta_t = far1_betas$beta1) {
  if (!is.matrix(X)) stop("far1_scan expects a matrix (grid x n)")
  L <- nrow(X); n <- ncol(X)
  if (n < 1L) {
    return(list(stat = 0, k = NA_integer_, m = NA_integer_,
                interval = c(NA_real_, NA_real_)))
  }
  s <- seq(0, 1, length.out = L)
  W <- (1 - beta_t(s)) * (X - rowMeans(X))
  CS <- cbind(0, t(apply(W, 1L, cumsum)))

  E_max <- 0; b_start <- NA_integer_; b_end <- NA_integer_
  for (d in 0:(n - 1L)) {
    m_idx <- seq.int(1L, n - d)
    k_idx <- m_idx + d
    block <- abs(CS[, k_idx + 1L, drop = FALSE] - CS[, m_idx, drop = FALSE])
    i <- which.max(block)
    val_d <- block[i] / (d + 1L)^gamma
    if (val_d > E_max) {
      E_max   <- val_d
      b_start <- ((i - 1L) %/% L) + 1L        # first column of the window
      b_end   <- b_start + d
    }
  }
  list(stat     = E_max * n^(-0.5 + gamma),
       k        = b_start - 1L,
       m        = b_end,
       interval = c(b_start - 1L, b_end) / n)
}

#' Beta operators from Birbilas (2025)
#'
#' The four pointwise multiplication operators \eqn{\beta(t)} of the
#' reference simulation study, as functions of `t`: `beta1` is
#' \eqn{0.7\sqrt{t}} (canonical), `beta2` its negative, `beta3` is
#' \eqn{0.7\sqrt{1-t}}, `beta4` its negative.
#'
#' @format A named list of four functions.
#' @examples
#' far1_betas$beta1(c(0, 0.5, 1))
#' @export
far1_betas <- list(
  beta1 = function(t)  0.7 * sqrt(t),
  beta2 = function(t) -0.7 * sqrt(t),
  beta3 = function(t)  0.7 * sqrt(1 - t),
  beta4 = function(t) -0.7 * sqrt(1 - t))

#' Null draws for the Birbilas statistic
#'
#' Monte-Carlo draws of the asymptotic limit \eqn{E_{q,\gamma}} of the
#' scaled statistic \eqn{n^{-1/2+\gamma} E_{n,\beta,\gamma}}, from the
#' shipped cache (1000 draws per gamma in \{0, 0.25, 0.4, 0.45\}).
#' Conservative for small `n`; for a finite-sample calibration simulate
#' the statistic under your H0 and pass the draws to [epidemic_test()]
#' via `null=`. Extend the asymptotic cache via
#' [simulate_far1_null()] + [save_null_draws()].
#'
#' @param gamma window-length exponent (scalar).
#' @return numeric vector of null draws.
#' @examples
#' quantile(far1_null_draws(gamma = 0.4), 0.95)
#' @export
far1_null_draws <- function(gamma) {
  stopifnot(length(gamma) == 1L)
  .null_get("far1", .null_key("far1", gamma = gamma),
            sprintf("simulate_far1_null(n_draws, gamma = %s)", gamma))
}

#' @rdname far1_null_draws
#' @param alpha significance level.
#' @param splits split-aggregation adjustment, see [critical_value()].
#' @return `far1_critical_value()`: scalar critical value.
#' @export
far1_critical_value <- function(gamma, alpha = 0.05, splits = 1) {
  critical_value(far1_null_draws(gamma), alpha, splits)
}
