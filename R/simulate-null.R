#' Simulate null draws for the scan statistics
#'
#' Monte-Carlo samplers of the four null distributions, for gammas,
#' kernels or baselines that the shipped caches do not cover. Store the
#' result with [save_null_draws()] to make the `*_null_draws()` accessors
#' - and [epidemic_test()] - pick it up automatically, or pass it
#' directly to `epidemic_test(null = )`.
#'
#' * `simulate_wilcoxon_null()`: the weighted Brownian-bridge range
#'   \eqn{\sup_{s<t} |B(t) - B(s)| \, w(|t-s|(1-|t-s|))} on a `grid`-point
#'   bridge. The shipped cache used `grid = 10000` and 30,000 draws per
#'   cell (weeks of compute); `grid = 2000` with a few thousand draws
#'   reproduces its quantiles to about 0.03.
#' * `simulate_lp_null()`: the weighted Lp functional of the centred
#'   Wiener sheet on a `grid` x `grid` lattice (shipped: `grid = 200`,
#'   300 draws).
#' * `simulate_rkhs_null()`: draws of [rkhs_statistic()] on i.i.d.
#'   samples of size `n` from the H0 sampler `rdist` (shipped:
#'   `n = 2000`, 1000 draws).
#' * `simulate_far1_null()`: the asymptotic limit
#'   \eqn{E_{q,\gamma}} on a `grid` x `grid` Brownian-motion lattice
#'   (shipped: `grid = 1200`, 1000 draws).
#'
#' For reproducibility set the RNG before calling: `set.seed()` with
#' `cores = 1`, or `RNGkind("L'Ecuyer-CMRG"); set.seed()` for parallel
#' runs. These are long computations at calibration scale; `progress`
#' prints a live done/total/ETA line.
#'
#' @param n_draws number of Monte-Carlo draws.
#' @param gamma weight exponent of the target cell.
#' @param beta Wilcoxon only: logarithmic weight exponent.
#' @param p Lp only: the Lp exponent.
#' @param kernel RKHS only: kernel name or function.
#' @param rdist RKHS only: H0 sampler, a function of `n` (e.g.
#'   [stats::rnorm]).
#' @param n RKHS only: sample size per draw.
#' @param grid discretisation resolution, see Details.
#' @param cores parallel workers (forked; use 1 on Windows).
#' @param progress print a progress line (default: only in interactive
#'   sessions).
#' @return numeric vector of `n_draws` null draws.
#' @seealso [save_null_draws()], [epidemic_test()]
#' @examples
#' set.seed(1)
#' # verification scale; calibration runs use far larger grid / n_draws
#' draws <- simulate_wilcoxon_null(20, gamma = 0.25, grid = 200)
#' quantile(draws, 0.95)
#' @export
simulate_wilcoxon_null <- function(n_draws, gamma, beta = 0, grid = 2000L,
                                   cores = 1L, progress = interactive()) {
  grid <- as.integer(grid)
  ls <- seq_len(grid)
  w_sym <- holder_weight(ls / grid * (1 - ls / grid), gamma, beta)
  one <- function(i) {
    W <- c(0, cumsum(rnorm(grid, 0, sqrt(1 / grid))))
    B <- W - (0:grid) / grid * W[grid + 1L]
    vals <- vapply(seq_len(grid), function(l) {
      max(abs(B[(l + 1L):(grid + 1L)] - B[seq_len(grid - l + 1L)])) * w_sym[l]
    }, numeric(1))
    max(vals)
  }
  unlist(.apply_progress(seq_len(n_draws), one, cores, progress))
}

#' @rdname simulate_wilcoxon_null
#' @export
simulate_lp_null <- function(n_draws, gamma, p = 2, grid = 200L, cores = 1L,
                             progress = interactive()) {
  grid <- as.integer(grid)
  N <- grid
  w <- holder_weight(seq_len(N - 1L) / N, gamma, 0)
  us <- seq_len(N) / N
  kl <- which(upper.tri(matrix(FALSE, N, N), diag = FALSE), arr.ind = TRUE)
  k_v <- kl[, 1L]; l_v <- kl[, 2L]
  lag_v <- l_v - k_v
  one <- function(i) {
    Z <- matrix(rnorm(N * N, 0, 1), N, N) / N
    W <- t(apply(apply(Z, 2L, cumsum), 1L, cumsum))   # 2-d prefix sums
    W_end_row <- W[N, ]
    W_end_col <- W[, N]
    W_end_end <- W[N, N]
    best <- 0
    for (p_idx in seq_along(k_v)) {
      k <- k_v[p_idx]; l <- l_v[p_idx]
      delta <- (W[, l] - W[, k]) - (us * W_end_row[l]) -
               ((l - k) / N) * W_end_col + (us * W_end_row[k]) +
               ((l - k) / N) * us * W_end_end
      int_u <- sum(abs(delta)^p) / N
      val   <- int_u^(1 / p) * w[lag_v[p_idx]]
      if (val > best) best <- val
    }
    best
  }
  unlist(.apply_progress(seq_len(n_draws), one, cores, progress))
}

#' @rdname simulate_wilcoxon_null
#' @export
simulate_rkhs_null <- function(n_draws, gamma, kernel = "gaussian",
                               rdist = stats::rnorm, n = 2000L, beta = 0,
                               cores = 1L, progress = interactive()) {
  kfn <- .resolve_kernel(kernel)
  one <- function(i) rkhs_statistic(rdist(n), gamma, kernel = kfn, beta = beta)
  unlist(.apply_progress(seq_len(n_draws), one, cores, progress))
}

#' @rdname simulate_wilcoxon_null
#' @export
simulate_far1_null <- function(n_draws, gamma, grid = 1200L, cores = 1L,
                                   progress = interactive()) {
  N <- as.integer(grid)
  one <- function(i) {
    x1 <- matrix(rnorm(N * N), nrow = N, ncol = N)
    x1 <- t(apply(x1, 1L, cumsum))                # N x N, each row a BM path
    CS <- rbind(0, apply(x1, 2L, cumsum))         # (N + 1) x N column prefixes
    E_max <- 0
    for (d in seq_len(N - 1L)) {
      m_idx <- seq.int(1L, N - d)
      k_idx <- m_idx + d
      block <- abs(CS[k_idx + 1L, , drop = FALSE] - CS[m_idx, , drop = FALSE])
      val_d <- max(block) / d^gamma
      if (val_d > E_max) E_max <- val_d
    }
    row_max <- max(abs(x1))                       # single-row boundary terms
    if (row_max > E_max) E_max <- row_max
    E_max * N^(-1 + gamma)
  }
  unlist(.apply_progress(seq_len(n_draws), one, cores, progress))
}
