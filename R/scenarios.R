# Pareto sampling and CDF via the closed-form inverse. Bit-identical to
# EnvStats::rpareto / EnvStats::ppareto (verified: qpareto inverts as
# location * (1 - p)^(-1/shape)), so the package needs no EnvStats
# dependency while reproducing the dissertation's draw streams exactly.
.rpareto <- function(n, location, shape) {
  location * (1 - runif(n))^(-1 / shape)
}

.ppareto <- function(q, location, shape) {
  ifelse(q >= location, 1 - (location / q)^shape, 0)
}

# Change-window indices for scalar scenarios. The dissertation's window
# I* = {k*+1, ..., k*+l} holds exactly l observations (exact_window =
# TRUE); the legacy code included both endpoints (l + 1 observations,
# exact_window = FALSE), kept for reproducing pre-audit runs.
.ec_index <- function(n, ec_interval, exact_window) {
  idx <- round(n * ec_interval)
  if (exact_window) (idx[1] + 1L):idx[2] else idx[1]:idx[2]
}

#' Benchmark scalar scenarios S1-S5
#'
#' The dissertation's scalar data-generating processes. Each draws an
#' i.i.d. sample of length `n` and, when `a != 0`, replaces the
#' observations in the change window with draws from the changed
#' distribution:
#'
#' * `generate_s1()`: N(0, 1) baseline; mean shift to N(a, 1).
#' * `generate_s2()`: squared-normal (chi-squared, df 1) baseline;
#'   in-window draws are `N(a, 1)^2`, i.e. noncentral with `ncp = a^2`.
#' * `generate_s3()`: Pareto(location 1, shape 3) baseline; location
#'   shifts to `1 + a`.
#' * `generate_s4()`: Beta(1, 9) baseline; shape shifts to
#'   Beta(`1 + a`, `9 - a`).
#' * `generate_s5()`: Pareto(1, 3) baseline; shape tightens to `3 + 10 a`.
#'
#' @param n sample length.
#' @param a change magnitude; `a = 0` is the null.
#' @param ec_interval two fractions in `[0, 1]`: endpoints of the change
#'   window as fractions of `n`.
#' @param exact_window if `TRUE` (default) the window is the
#'   dissertation's \eqn{I^* = \{k^*+1, \dots, k^*+l\}} with exactly
#'   `l` observations; `FALSE` reproduces the legacy convention that
#'   included both endpoints (`l + 1` observations).
#' @return numeric vector of length `n`.
#' @seealso [scenario_cdf()] for the matching null CDFs,
#'   [generate_f1()] for the functional scenarios.
#' @examples
#' set.seed(1)
#' x <- generate_s1(200, a = 1)
#' epidemic_test(x, "wilcoxon", gamma = 0.25)
#' @export
generate_s1 <- function(n, a = 0, ec_interval = c(0.4, 0.7),
                        exact_window = TRUE) {
  X <- rnorm(n, 0, 1)
  if (a != 0) {
    k <- .ec_index(n, ec_interval, exact_window)
    X[k] <- rnorm(length(k), a, 1)
  }
  X
}

#' @rdname generate_s1
#' @export
generate_s2 <- function(n, a = 0, ec_interval = c(0.4, 0.7),
                        exact_window = TRUE) {
  X <- rnorm(n, 0, 1)^2
  if (a != 0) {
    k <- .ec_index(n, ec_interval, exact_window)
    X[k] <- rnorm(length(k), a, 1)^2
  }
  X
}

#' @rdname generate_s1
#' @export
generate_s3 <- function(n, a = 0, ec_interval = c(0.4, 0.7),
                        exact_window = TRUE) {
  X <- .rpareto(n, location = 1, shape = 3)
  if (a != 0) {
    k <- .ec_index(n, ec_interval, exact_window)
    X[k] <- .rpareto(length(k), location = 1 + a, shape = 3)
  }
  X
}

#' @rdname generate_s1
#' @export
generate_s4 <- function(n, a = 0, ec_interval = c(0.4, 0.7),
                        exact_window = TRUE) {
  X <- rbeta(n, 1, 9)
  if (a != 0) {
    k <- .ec_index(n, ec_interval, exact_window)
    X[k] <- rbeta(length(k), 1 + a, 9 - a)
  }
  X
}

#' @rdname generate_s1
#' @export
generate_s5 <- function(n, a = 0, ec_interval = c(0.4, 0.7),
                        exact_window = TRUE) {
  X <- .rpareto(n, location = 1, shape = 3)
  if (a != 0) {
    k <- .ec_index(n, ec_interval, exact_window)
    X[k] <- .rpareto(length(k), location = 1, shape = 3 + 10 * a)
  }
  X
}

#' Benchmark functional scenarios F1-F5
#'
#' The dissertation's functional data-generating processes. Each returns a
#' matrix with `l + 1` rows (the evaluation grid on `[0, 1]`) and `n`
#' columns (one curve per observation):
#'
#' * `generate_f1()`: cosine-basis curves
#'   \eqn{X_j(t) = \sqrt{2} \sum_{k=1}^{10} \xi_{jk} \cos(k\pi t)/k} with
#'   N(0, 1) coefficients under H0 and centred-Pareto differences
#'   (`Pa(1, a) - Pa(1, a)`) in the change window.
#' * `generate_f2()`: same basis, Pareto-difference coefficients
#'   everywhere; the in-window Pareto shape is `a`, elsewhere 1.5, so
#'   `a = 1.5` is a no-change configuration and `a = 0` the conventional
#'   null.
#' * `generate_f3()`: integrated partial-sum (random polygonal) curves
#'   with the same coefficient change as F1.
#' * `generate_f4()`: the FAR(1) functional time series of Birbilas
#'   (2025), \eqn{X_i = \beta(t) X_{i-1} + g(i/n) + Y_i} with Wiener
#'   innovations and epidemic jump \eqn{g = a\sqrt{t}} on the window;
#'   note the observations are serially dependent, and the window
#'   convention is exact by construction (no `exact_window` argument).
#' * `generate_f5()`: i.i.d. standard Wiener (random-walk) curves; in the
#'   change window each curve receives an additive **spike**
#'   \eqn{a \, f_\alpha(t)} with the dissertation's spike function
#'   \deqn{f_\alpha(t) = \big(1 - |t - 1/2|^\alpha\big)^{1/\alpha},
#'     \qquad \alpha \in (0, 1/2),}
#'   which peaks at \eqn{f_\alpha(1/2) = 1} with a Hoelder-\eqn{\alpha}
#'   cusp (the function of the Faber-Schauder smoothing figure; default
#'   `spike_alpha = 0.45` as printed there). Restored from the paper-1
#'   draft's scenario (S2) (archive ids `p1.TD.S2.spike` /
#'   `wiener_TD2_spike`), which was cut before publication; the archived
#'   forms used a tent-power bump \eqn{(1 - |2t - 1|)^k} (\eqn{k = 4} in
#'   the draft text, \eqn{k = 20} in the code), replaced here by the
#'   dissertation's \eqn{f_\alpha} parametrisation. The draft smoothed
#'   the curves with the Faber-Schauder basis before testing; pair with
#'   [reconstruct()] to replicate that pipeline. This is the scenario
#'   behind the claim that point-evaluation projections (e.g.
#'   [project_c_max_var()]) are particularly informative for spike-type
#'   (time-localised) changes.
#'
#' @inheritParams generate_s1
#' @param l grid resolution: curves are sampled at `l + 1` equispaced
#'   points on `[0, 1]`.
#' @param beta_t F4 only: the FAR(1) operator \eqn{\beta(t)}, see
#'   [far1_betas].
#' @param spike_alpha F5 only: the Hoelder exponent \eqn{\alpha} of the
#'   spike function \eqn{f_\alpha}; the dissertation uses
#'   \eqn{\alpha \in (0, 1/2)}, smaller values giving a sharper cusp.
#' @return numeric matrix, `(l + 1) x n`.
#' @seealso [project_l2()] and friends to reduce curves to scalars,
#'   [far1_statistic()] which consumes the matrix directly.
#' @examples
#' set.seed(1)
#' X <- generate_f1(80, a = 3)
#' epidemic_test(project_l2(X), "rkhs", gamma = 0.25)
#' @export
generate_f1 <- function(n, a = 0, ec_interval = c(0.4, 0.7), l = 100,
                        exact_window = TRUE) {
  s <- (0:l) / l
  X <- matrix(0, nrow = l + 1L, ncol = n)
  idx <- round(n * ec_interval)
  j_lo <- if (exact_window) idx[1] + 1L else idx[1]
  in_window <- function(j) a != 0 && j >= j_lo && j <= idx[2]
  for (j in seq_len(n)) {
    X[, j] <- .f1_event(s, window = in_window(j), a = a)
  }
  X
}

.f1_event <- function(t.s, window, a) {
  d <- 10L
  coefs <- if (window) {
    .rpareto(d, location = 1, shape = a) - .rpareto(d, location = 1, shape = a)
  } else {
    rnorm(d, 0, 1)
  }
  out <- numeric(length(t.s))
  for (k in seq_len(d)) out <- out + coefs[k] * cos(k * pi * t.s) / k
  out * sqrt(2)
}

#' @rdname generate_f1
#' @export
generate_f2 <- function(n, a = 0, ec_interval = c(0.4, 0.7), l = 100,
                        exact_window = TRUE) {
  s <- (0:l) / l
  X <- matrix(0, nrow = l + 1L, ncol = n)
  idx <- round(n * ec_interval)
  j_lo <- if (exact_window) idx[1] + 1L else idx[1]
  in_window <- function(j) a != 0 && j >= j_lo && j <= idx[2]
  for (j in seq_len(n)) {
    theta <- if (in_window(j)) a else 1.5
    X[, j] <- .f2_event(s, theta)
  }
  X
}

.f2_event <- function(t.s, theta) {
  d <- 10L
  eta <- 1
  coefs <- .rpareto(d, location = eta, shape = theta) -
           .rpareto(d, location = eta, shape = theta)
  out <- numeric(length(t.s))
  for (k in seq_len(d)) out <- out + coefs[k] * cos(k * pi * t.s) / k
  out * sqrt(2)
}

#' @rdname generate_f1
#' @export
generate_f3 <- function(n, a = 0, ec_interval = c(0.4, 0.7), l = 100,
                        exact_window = TRUE) {
  s <- (0:l) / l
  X <- matrix(0, nrow = l + 1L, ncol = n)
  idx <- round(n * ec_interval)
  j_lo <- if (exact_window) idx[1] + 1L else idx[1]
  in_window <- function(j) a != 0 && j >= j_lo && j <= idx[2]
  for (j in seq_len(n)) {
    X[, j] <- .f3_event(s, window = in_window(j), a = a)
  }
  X
}

.f3_event <- function(t.s, window, a) {
  d <- 10L
  xi <- if (window) {
    .rpareto(d, location = 1, shape = a) - .rpareto(d, location = 1, shape = a)
  } else {
    rnorm(d, 0, 1)
  }
  cs <- c(0, cumsum(xi))                       # cs[k+1] = xi_1 + ... + xi_k
  dt <- d * t.s
  k_floor <- pmin(as.integer(floor(dt)), d)    # cap at d for t == 1
  frac    <- dt - k_floor
  part <- cs[k_floor + 1L]
  extra_idx <- pmin(k_floor + 1L, d)
  extra <- frac * xi[extra_idx]
  part + extra
}

#' @rdname generate_f1
#' @export
generate_f4 <- function(n, a = 0, ec_interval = c(0.25, 0.45), l = 100,
                        beta_t = far1_betas$beta1) {
  s <- seq(0, 1, length.out = l + 1L)
  beta_grid <- beta_t(s)
  delta_t <- 1 / l

  k_lo <- round(n * ec_interval[1])
  k_hi <- round(n * ec_interval[2])

  X <- matrix(0, nrow = l + 1L, ncol = n)
  X_prev <- numeric(l + 1L)            # X_0 = 0
  for (i in seq_len(n)) {
    Y_i <- c(0, cumsum(sqrt(delta_t) * rnorm(l, 0, 1)))
    X_i <- beta_grid * X_prev + Y_i
    if (a != 0 && i > k_lo && i <= k_hi) {
      X_i <- X_i + a * sqrt(s)
    }
    X[, i] <- X_i
    X_prev <- X_i
  }
  X
}

#' @rdname generate_f1
#' @export
generate_f5 <- function(n, a = 0, ec_interval = c(0.4, 0.7), l = 100,
                        spike_alpha = 0.45, exact_window = TRUE) {
  s <- (0:l) / l
  spike <- (1 - abs(s - 0.5)^spike_alpha)^(1 / spike_alpha)
  delta_t <- 1 / l
  X <- matrix(0, nrow = l + 1L, ncol = n)
  idx <- round(n * ec_interval)
  j_lo <- if (exact_window) idx[1] + 1L else idx[1]
  in_window <- function(j) a != 0 && j >= j_lo && j <= idx[2]
  for (j in seq_len(n)) {
    W <- c(0, cumsum(sqrt(delta_t) * rnorm(l, 0, 1)))
    X[, j] <- if (in_window(j)) W + a * spike else W
  }
  X
}

#' Null CDF of a benchmark scenario
#'
#' Returns the pointwise null CDF `F_0` matching a scenario, for use as
#' the `cdf` argument of [lp_statistic()] / [epidemic_test()]: standard
#' normal (S1), chi-squared df 1 (S2), Pareto(1, 3) (S3, S5), Beta(1, 9)
#' (S4). The functional scenarios (F1-F5) have no closed-form CDF after
#' projection; for them the returned function applies the empirical CDF
#' of the sample itself, which is also what `cdf = NULL` does.
#'
#' @param scenario scenario id: `"S1"` ... `"S5"`, `"F1"` ... `"F5"`
#'   (case-insensitive).
#' @return a vectorised function `F_0(x)`.
#' @examples
#' set.seed(1)
#' x <- generate_s1(150, a = 1)
#' lp_statistic(x, gamma = 0.25, cdf = scenario_cdf("S1"))
#' @export
scenario_cdf <- function(scenario) {
  switch(toupper(as.character(scenario)),
    S1 = function(x) pnorm(x, 0, 1),
    S2 = function(x) pchisq(x, df = 1, ncp = 0),
    S3 = ,
    S5 = function(x) .ppareto(x, location = 1, shape = 3),
    S4 = function(x) pbeta(x, shape1 = 1, shape2 = 9),
    F1 = ,
    F2 = ,
    F3 = ,
    F4 = ,
    F5 = function(x) ecdf(x)(x),
    stop("unknown scenario: ", scenario))
}
