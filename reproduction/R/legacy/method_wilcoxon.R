#' Antisymmetric-kernel (Wilcoxon) scan test (Chapter 3 of the dissertation).
#'
#' For each window [k, m] computes
#'   H(k,m) = sum_{i in [k,m], j not in [k,m]} sign(f(X_i) - f(X_j))
#' where f is the feature map selected by `kernel`:
#'   - "h_w" (default): f = identity -- raw Wilcoxon kernel h^w(x,y) =
#'     sign(x - y). Sensitive to median / mean shifts.
#'   - "h_0"          : f = abs      -- norm-based Wilcoxon kernel
#'     h_0(x,y) = sign(|x| - |y|). Sensitive to tail / scale changes;
#'     this is what paper 1 (A1) uses for its functional simulations
#'     (see `helper_functions.R::get.h.Wilcoxon`, line 571).
#'
#' Reports the symmetrically-weighted scan statistic
#'   T = max_{k < m} |H(k,m)| / w(l/n)    with l = m - k
#'       / ( n^{1.5} * sigma_hat )
#' using the symmetric weight w(l/n) = rho(l(n-l)/n^2).
#' sigma_hat^2 = mean_j ( mean_i sign(...) )^2 -- feature-invariant for
#' continuous symmetric distributions, so the same null distribution
#' (`.raw/T2_gamma_beta.csv`) applies to both kernel choices.
#'
#' Callers load shared helpers once (integral_image.R, weights.R, null_cache.R).

wilcoxon_statistic <- function(X, gamma, beta = 0, kernel = "h_w", ...) {
  feature <- switch(kernel,
    h_w = identity,
    h_0 = abs,
    stop("unknown wilcoxon kernel: ", kernel))
  X <- feature(X)
  n <- length(X)
  H <- outer(X, X, FUN = function(x, y) sign(x - y))
  sigma <- sqrt(mean(colMeans(H)^2))
  if (!is.finite(sigma) || sigma == 0) return(0)
  w <- weight_cache(n, gamma, beta, symmetric = TRUE)

  # Vectorised scan: enumerate all (k, m) with k < m, then compute
  # `sumRegion(sat, k, 1, m, n) - sumRegion(sat, k, k, m, m)` over all pairs
  # via matrix-index lookups on the integral image. Matches the previous
  # nested loop bit-exactly; ~17x faster at n = 300.
  sat   <- NumMatrix(H)$sat
  max_r <- nrow(sat); max_c <- ncol(sat)
  km    <- which(upper.tri(matrix(FALSE, n, n), diag = FALSE), arr.ind = TRUE)
  k_vec <- km[, 1L]; m_vec <- km[, 2L]
  r2    <- pmin(m_vec + 2L, max_r)
  c2    <- pmin(m_vec + 2L, max_c)
  k1    <- k_vec + 1L
  c_end <- rep(max_c, length(k_vec))
  two   <- rep(2L,    length(k_vec))
  R_km  <- sat[cbind(r2, c_end)] - sat[cbind(r2, two)] -
           sat[cbind(k1, c_end)] + sat[cbind(k1, two)]
  B_km  <- sat[cbind(r2, c2)]    - sat[cbind(r2, k1)]  -
           sat[cbind(k1, c2)]    + sat[cbind(k1, k1)]
  vals  <- abs(R_km - B_km) / w[m_vec - k_vec]
  vals[!is.finite(vals)] <- -Inf
  max(vals) / (n^1.5 * sigma)
}

#' Corrected scan: the statistic exactly as defined in the header above.
#'
#' `wilcoxon_statistic` / `wilcoxon_scan` carry two off-by-one errors inherited
#' from mixing `sumRegion`'s 0-based inclusive indexing with 1-based matrix
#' indices. The row window is evaluated at rows (k+1)..(m+1) while being
#' reported as (k, m), which is a harmless relabelling, and column 1 is dropped
#' from the complement, which is not: every H(i, 1) term is missing, so the
#' statistic is systematically understated (about 0.6 to 2.3 per cent for
#' n between 50 and 400, verified against brute force). Because the critical
#' values come from a cached asymptotic null rather than from the same code
#' path, the bias is not cancelled by calibration; it makes the test slightly
#' conservative.
#'
#' The original archive implementation (`helper_functions.R:598`,
#' `sum(h.cache[k:m, ]) - sum(h.cache[k:m, k:m])`) is correct, and this
#' function reproduces it. It is provided alongside rather than in place of the
#' existing functions so that every previously published simulation result
#' stays bit-reproducible; callers opt in.
#' Gamma-independent core shared by the three exact entry points. Returns the
#' unweighted |H(k,m)| over every window together with the window lengths, the
#' sigma normaliser and the (k, m) index vectors.
.wilcoxon_core_exact <- function(X, kernel) {
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

  km    <- which(upper.tri(matrix(FALSE, n, n), diag = FALSE), arr.ind = TRUE)
  k_vec <- km[, 1L]; m_vec <- km[, 2L]
  rows  <- P[cbind(m_vec + 1L, n + 1L)] - P[cbind(k_vec, n + 1L)]
  block <- P[cbind(m_vec + 1L, m_vec + 1L)] - P[cbind(m_vec + 1L, k_vec)] -
           P[cbind(k_vec,      m_vec + 1L)] + P[cbind(k_vec,      k_vec)]
  list(diff_abs = abs(rows - block), l_vec = m_vec - k_vec,
       k_vec = k_vec, m_vec = m_vec, sigma = sigma, n = n)
}

wilcoxon_scan_exact <- function(X, gamma, beta = 0, kernel = "h_w", ...) {
  co <- .wilcoxon_core_exact(X, kernel)
  if (is.null(co)) {
    return(list(stat = 0, k = NA_integer_, m = NA_integer_,
                interval = c(NA_real_, NA_real_)))
  }
  w    <- weight_cache(co$n, gamma, beta, symmetric = TRUE)
  vals <- co$diff_abs / w[co$l_vec]
  vals[!is.finite(vals)] <- -Inf
  i_max <- which.max(vals)
  list(
    stat     = vals[i_max] / (co$n^1.5 * co$sigma),
    k        = co$k_vec[i_max],
    m        = co$m_vec[i_max],
    interval = c(co$k_vec[i_max], co$m_vec[i_max]) / co$n
  )
}

#' Exact counterparts of `wilcoxon_statistic` and `wilcoxon_statistic_gammas`,
#' for power runs that opt out of the legacy off-by-one. `driver.R` selects
#' between the two families through its `WILCOXON_SCAN` constant.
wilcoxon_statistic_exact <- function(X, gamma, beta = 0, kernel = "h_w", ...) {
  co <- .wilcoxon_core_exact(X, kernel)
  if (is.null(co)) return(0)
  w    <- weight_cache(co$n, gamma, beta, symmetric = TRUE)
  vals <- co$diff_abs / w[co$l_vec]
  vals[!is.finite(vals)] <- -Inf
  max(vals) / (co$n^1.5 * co$sigma)
}

wilcoxon_statistic_gammas_exact <- function(X, gammas, beta = 0,
                                            kernel = "h_w", ...) {
  co <- .wilcoxon_core_exact(X, kernel)
  if (is.null(co)) return(rep(0, length(gammas)))
  norm <- co$n^1.5 * co$sigma
  vapply(gammas, function(g) {
    w    <- weight_cache(co$n, g, beta, symmetric = TRUE)
    vals <- co$diff_abs / w[co$l_vec]
    vals[!is.finite(vals)] <- -Inf
    max(vals) / norm
  }, numeric(1))
}

#' CORRECTED scan family: the statistic as the dissertation defines it.
#'
#' Beyond the summed-area off-by-one that `wilcoxon_scan_exact` fixed, the
#' legacy family carries two deeper defects, found in the 2026-09 audit:
#'
#'   1. WEIGHT DIRECTION. `weight_cache` returns w = 1/rho (its docstring says
#'      so). The thesis statistic DIVIDES |Delta| by rho, i.e. MULTIPLIES by
#'      w; the Lp and RKHS statistics and the Wilcoxon NULL GENERATOR
#'      (`generators/wilcoxon_null.R:71`, and the GPU script that built
#'      `.raw/T2_gamma_beta.csv`) all multiply. Every legacy Wilcoxon entry
#'      point divides, so at gamma > 0 the statistic shrinks while its
#'      critical value grows: under H_0 at n = 400 the empirical 95% point is
#'      1.17 against a cached CV of 2.54 at gamma = 0.25, and rejection is
#'      arithmetically near-impossible. This -- not conservative critical
#'      values -- is why the printed tables report zero Wilcoxon power at
#'      gamma > 0. Inherited from `archive/.../helper_functions.R:599` (the
#'      CUSUM in the same file multiplies at :553).
#'   2. WEIGHT INDEX. The scanned window holds l = b - a + 1 observations but
#'      the weight is looked up at l - 1; the null pairs a lag-l bridge
#'      increment with weight at l/n. Off by one, worst at short windows
#'      (factor 2^gamma at l = 2).
#'
#' The corrected family multiplies by the weight, indexes it by the true
#' window length, and scans every window 1 <= a <= b <= n including length-1
#' windows (the thesis allows m = k + 1). Normalisation and sigma_hat follow
#' Theorem main-teo2 unchanged. Legacy and exact families are kept verbatim
#' so published numbers stay reproducible; select via WILCOXON_SCAN.
.wilcoxon_core_corrected <- function(X, kernel) {
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

  # Windows I = {a..b} for all 1 <= a <= b <= n. The full window (1, n) has an
  # empty complement, so Delta = 0 and w[n] = 0; it drops out on its own.
  b_vec <- rep.int(seq_len(n), seq_len(n))
  a_vec <- sequence(seq_len(n))
  rows  <- P[cbind(b_vec + 1L, n + 1L)] - P[cbind(a_vec, n + 1L)]
  block <- P[cbind(b_vec + 1L, b_vec + 1L)] - P[cbind(b_vec + 1L, a_vec)] -
           P[cbind(a_vec, b_vec + 1L)] + P[cbind(a_vec, a_vec)]
  list(diff_abs = abs(rows - block), l_vec = b_vec - a_vec + 1L,
       a_vec = a_vec, b_vec = b_vec, sigma = sigma, n = n)
}

wilcoxon_statistic_corrected <- function(X, gamma, beta = 0, kernel = "h_w", ...) {
  co <- .wilcoxon_core_corrected(X, kernel)
  if (is.null(co)) return(0)
  w    <- weight_cache(co$n, gamma, beta, symmetric = TRUE)
  vals <- co$diff_abs * w[co$l_vec]
  max(vals) / (co$n^1.5 * co$sigma)
}

wilcoxon_statistic_gammas_corrected <- function(X, gammas, beta = 0,
                                                kernel = "h_w", ...) {
  co <- .wilcoxon_core_corrected(X, kernel)
  if (is.null(co)) return(rep(0, length(gammas)))
  norm <- co$n^1.5 * co$sigma
  vapply(gammas, function(g) {
    w <- weight_cache(co$n, g, beta, symmetric = TRUE)
    max(co$diff_abs * w[co$l_vec]) / norm
  }, numeric(1))
}

wilcoxon_scan_corrected <- function(X, gamma, beta = 0, kernel = "h_w", ...) {
  co <- .wilcoxon_core_corrected(X, kernel)
  if (is.null(co)) {
    return(list(stat = 0, k = NA_integer_, m = NA_integer_,
                interval = c(NA_real_, NA_real_)))
  }
  w     <- weight_cache(co$n, gamma, beta, symmetric = TRUE)
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

#' Null distribution: precomputed 2-sided Hoelder-type draws from
#' `.raw/T2_gamma_beta.csv`. Indexed by (gamma, beta); `beta` only varies at
#' gamma = 0.5 and is forced to 0 elsewhere. `exact_gamma0` swaps the
#' (gamma = 0, beta = 0) cell for the closed-form Kuiper law -- the cached
#' column for that cell is provably wrong (its 95% point, 1.770, exceeds the
#' continuous-limit value 1.7473, and its maximum over 30,000 draws sits below
#' the exact 99.9% point); printed Table 3.1 matches the closed form.
wilcoxon_null_distribution <- function(n, gamma, beta = 0,
                                       exact_gamma0 = FALSE, ...) {
  load_wilcoxon_null(gamma, beta, exact_gamma0 = exact_gamma0)
}

#' Critical value for the Wilcoxon scan at level alpha (with m.d-correction).
#' Wraps `wilcoxon_null_distribution` + `critical_value` in one call so callers
#' don't have to know about the precomputed-null pipeline. `n` is forwarded
#' for signature parity with other methods even though the cached asymptotic
#' null does not depend on it.
wilcoxon_critical_value <- function(n, gamma, beta = 0, alpha = 0.05, m.d = 1L,
                                    exact_gamma0 = FALSE) {
  critical_value(wilcoxon_null_distribution(n, gamma, beta,
                                            exact_gamma0 = exact_gamma0),
                 alpha, m.d)
}

#' Wilcoxon scan statistic at multiple gammas, sharing the gamma-independent
#' work. About 70% of `wilcoxon_statistic` (the H matrix, integral image,
#' R_km/B_km lookups, sigma) does not depend on gamma; calling
#' `wilcoxon_statistic` once per gamma in a hot loop redoes that work for
#' every gamma. This batched variant computes the gamma-independent piece
#' once and then folds in the gamma-specific weight + max for each gamma.
#'
#' Returns numeric(length(gammas)) -- the per-gamma scan statistics, in the
#' same order as `gammas`. Numerically equivalent to
#' `vapply(gammas, function(g) wilcoxon_statistic(X, g, beta, kernel), 0)`.
wilcoxon_statistic_gammas <- function(X, gammas, beta = 0, kernel = "h_w", ...) {
  feature <- switch(kernel,
    h_w = identity,
    h_0 = abs,
    stop("unknown wilcoxon kernel: ", kernel))
  X <- feature(X)
  n <- length(X)
  H <- outer(X, X, FUN = function(x, y) sign(x - y))
  sigma <- sqrt(mean(colMeans(H)^2))
  if (!is.finite(sigma) || sigma == 0) return(rep(0, length(gammas)))

  sat   <- NumMatrix(H)$sat
  max_r <- nrow(sat); max_c <- ncol(sat)
  km    <- which(upper.tri(matrix(FALSE, n, n), diag = FALSE), arr.ind = TRUE)
  k_vec <- km[, 1L]; m_vec <- km[, 2L]
  r2    <- pmin(m_vec + 2L, max_r)
  c2    <- pmin(m_vec + 2L, max_c)
  k1    <- k_vec + 1L
  c_end <- rep(max_c, length(k_vec))
  two   <- rep(2L,    length(k_vec))
  R_km  <- sat[cbind(r2, c_end)] - sat[cbind(r2, two)] -
           sat[cbind(k1, c_end)] + sat[cbind(k1, two)]
  B_km  <- sat[cbind(r2, c2)]    - sat[cbind(r2, k1)]  -
           sat[cbind(k1, c2)]    + sat[cbind(k1, k1)]
  diff_abs <- abs(R_km - B_km)
  l_vec    <- m_vec - k_vec
  norm     <- n^1.5 * sigma

  vapply(gammas, function(g) {
    w    <- weight_cache(n, g, beta, symmetric = TRUE)
    vals <- diff_abs / w[l_vec]
    vals[!is.finite(vals)] <- -Inf
    max(vals) / norm
  }, numeric(1))
}

#' Full Wilcoxon scan: max statistic plus argmax and detected interval.
#'
#' Same computation as `wilcoxon_statistic`; additionally returns the
#' (k, m) pair that achieves the max and the corresponding fractional
#' interval (k/n, m/n). Use this when the detected change window is
#' the quantity of interest (e.g. the Litgrid application); for power
#' simulations prefer the lighter `wilcoxon_statistic`.
wilcoxon_scan <- function(X, gamma, beta = 0, kernel = "h_w", ...) {
  feature <- switch(kernel,
    h_w = identity, h_0 = abs,
    stop("unknown wilcoxon kernel: ", kernel))
  X <- feature(X)
  n <- length(X)
  H <- outer(X, X, FUN = function(x, y) sign(x - y))
  sigma <- sqrt(mean(colMeans(H)^2))
  if (!is.finite(sigma) || sigma == 0) {
    return(list(stat = 0, k = NA_integer_, m = NA_integer_,
                interval = c(NA_real_, NA_real_)))
  }
  w <- weight_cache(n, gamma, beta, symmetric = TRUE)
  sat   <- NumMatrix(H)$sat
  max_r <- nrow(sat); max_c <- ncol(sat)
  km    <- which(upper.tri(matrix(FALSE, n, n), diag = FALSE), arr.ind = TRUE)
  k_vec <- km[, 1L]; m_vec <- km[, 2L]
  r2    <- pmin(m_vec + 2L, max_r)
  c2    <- pmin(m_vec + 2L, max_c)
  k1    <- k_vec + 1L
  c_end <- rep(max_c, length(k_vec))
  two   <- rep(2L,    length(k_vec))
  R_km  <- sat[cbind(r2, c_end)] - sat[cbind(r2, two)] -
           sat[cbind(k1, c_end)] + sat[cbind(k1, two)]
  B_km  <- sat[cbind(r2, c2)]    - sat[cbind(r2, k1)]  -
           sat[cbind(k1, c2)]    + sat[cbind(k1, k1)]
  vals  <- abs(R_km - B_km) / w[m_vec - k_vec]
  vals[!is.finite(vals)] <- -Inf
  i_max <- which.max(vals)
  list(
    stat     = vals[i_max] / (n^1.5 * sigma),
    k        = k_vec[i_max],
    m        = m_vec[i_max],
    interval = c(k_vec[i_max], m_vec[i_max]) / n
  )
}
