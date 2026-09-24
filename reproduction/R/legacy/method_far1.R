#' Birbilas's epidemic-change test E_{n,beta,gamma}.
#'
#' Reference: dissertation eq. `test_stat_2`
#' (`archive/phd/example_algimantas_far1_dissertation.tex`):
#'
#'   E_{n,beta,gamma} = max_{0 <= k < m <= n} max_{t in [0,1]}
#'                        |(I - beta(t)) sum_{i = k+1}^m (X_i(t) - bar_X_n(t))|
#'                        / (m - k)^gamma
#'
#'   bar_X_n(t) = n^{-1} sum_i X_i(t)
#'
#' Reject H_0: g = 0 at level alpha when n^{-1/2 + gamma} E > C''_alpha,
#' where C''_alpha is the (1 - alpha) quantile of the asymptotic limit
#' E_{q,gamma} (cached in `.raw/E_q_gamma.csv` via `far1_null.R`).
#'
#' Birbilas's defaults: beta_1(t) = 0.7 sqrt(t), gamma = 0.4. He uses
#' pointwise beta operators throughout the simulation study.
#'
#' Input contract: numeric matrix with l + 1 rows (grid points on [0, 1])
#' and n columns (one functional observation per column). Wire via
#' `make_run(matrix_in = TRUE)` so smoke does not fall through the
#' per-row scalar path.

far1_statistic <- function(X, gamma,
                               beta_t = function(t) 0.7 * sqrt(t), ...) {
  if (!is.matrix(X)) stop("far1_statistic expects matrix(l+1, n)")
  L <- nrow(X); n <- ncol(X)
  if (n < 2L) return(0)
  s <- seq(0, 1, length.out = L)
  W <- (1 - beta_t(s)) * (X - rowMeans(X))   # (I - beta) on centered curves

  # Column prefix sums per row: CS[, j + 1] = sum_{i=1..j} W[, i]. Built via
  # row-cumsum of W (apply with MARGIN=1 returns n x L; transpose to L x n)
  # then prepend a zero column.
  CS <- cbind(0, t(apply(W, 1L, cumsum)))    # L x (n + 1)

  # Iterate over lag d = k - m. For fixed d, the Birbilas-notation pairs
  # are m in 1..n-d, k = m + d. The sum from m to k inclusive is (d + 1)
  # columns, so block = CS[, k + 1] - CS[, m] (using CS[, 1] = 0 boundary).
  # Per-d block size is L x (n - d) -- peak L x (n - 1) at d = 1.
  E_max <- 0
  for (d in seq_len(n - 1L)) {
    m_idx <- seq.int(1L, n - d)
    k_idx <- m_idx + d
    block <- abs(CS[, k_idx + 1L, drop = FALSE] - CS[, m_idx, drop = FALSE])
    val_d <- max(block) / d^gamma
    if (val_d > E_max) E_max <- val_d
  }
  E_max * n^(-0.5 + gamma)
}

#' CORRECTED statistic (2026-09 audit, B1/B2): the printed denominator.
#'
#' `far1_statistic` sums a window of d + 1 columns but divides by
#' d^gamma, i.e. (length - 1)^gamma where the printed eq:far1-stat divides
#' the (m - k)-term sum by (m - k)^gamma -- +21% at Birbilas's canonical
#' gamma = 0.4 -- and it never scans length-1 windows (the print allows
#' m = k + 1). Self-cancelling for LEVEL under the default adjusted CV, whose
#' finite null is simulated with the same function, so the printed F4 rows
#' are a validly calibrated but different test; under the asymptotic CV the
#' legacy version over-rejects. The corrected variant pairs the regenerated
#' finite null (`load_far1_finite_null(..., variant = "corrected")`).
far1_statistic_corrected <- function(X, gamma,
                                         beta_t = function(t) 0.7 * sqrt(t),
                                         ...) {
  if (!is.matrix(X)) stop("far1_statistic_corrected expects matrix(l+1, n)")
  L <- nrow(X); n <- ncol(X)
  if (n < 1L) return(0)
  s <- seq(0, 1, length.out = L)
  W <- (1 - beta_t(s)) * (X - rowMeans(X))
  CS <- cbind(0, t(apply(W, 1L, cumsum)))

  E_max <- 0
  for (d in 0:(n - 1L)) {                     # window of d + 1 columns
    m_idx <- seq.int(1L, n - d)
    k_idx <- m_idx + d
    block <- abs(CS[, k_idx + 1L, drop = FALSE] - CS[, m_idx, drop = FALSE])
    val_d <- max(block) / (d + 1L)^gamma      # printed: length^gamma
    if (val_d > E_max) E_max <- val_d
  }
  E_max * n^(-0.5 + gamma)
}

#' Null distribution of n^{-1/2 + gamma} E_{n,beta,gamma}: precomputed
#' Wiener-sheet draws from `.raw/E_q_gamma.csv` (generator:
#' `generators/far1_null.R`). Indexed by gamma. Conservative for
#' n <= 500 -- prefer `far1_finite_null` for those.
far1_null_distribution <- function(n, gamma, ...) {
  load_far1_null(gamma)
}

#' Adjusted (finite-sample) null: empirical (1 - alpha)-quantile of
#' n^{-1/2 + gamma} E_{n,beta,gamma} simulated under H_0 at the target n.
#' Reads from `.raw/E_n_beta_gamma_finite.csv` (generator:
#' `generators/far1_finite_null.R`).
#'
#' NOTE (audit 2026-09, B4): this is an ORACLE quantile -- H_0 samples drawn
#' from the true generator with the true beta -- not Birbilas's Algorithm
#' `alg:zero`, which block-bootstraps the OBSERVED sample with an estimated
#' beta-hat and adjusts the level per eq:adjusted_fit3. The two agree in the
#' limit but are different procedures, and the oracle version is unavailable
#' to a practitioner. It is also asymmetric across the F4 comparison: Birbilas
#' and RKHS get finite-sample F4 quantiles while Lp and Wilcoxon read
#' asymptotic caches.
far1_finite_null <- function(n, gamma, beta_id = "beta1",
                                 variant = "legacy", ...) {
  load_far1_finite_null(n, gamma, beta_id, variant = variant)
}

#' Registry of beta operators from dissertation eq:my_betas. Use a string
#' id (beta1..beta4) so the finite-null cache can key on it. beta1 is
#' canonical for power simulations.
far1_beta_registry <- list(
  beta1 = function(t)  0.7 * sqrt(t),
  beta2 = function(t) -0.7 * sqrt(t),
  beta3 = function(t)  0.7 * sqrt(1 - t),
  beta4 = function(t) -0.7 * sqrt(1 - t))
