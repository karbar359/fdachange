# Brute-force gate: each statistic must equal a literal implementation of
# its printed formula to machine precision (ported from the thesis audit's
# verification gate).

TOL <- 1e-10
SEED <- 7L

test_that("wilcoxon_statistic equals the printed formula, brute-forced", {
  # Printed: windows I(k, m) = {k+1..m}, l = m - k in 1..n, MULTIPLY the
  # increment by the symmetric weight indexed at the true window length.
  brute_wilcoxon <- function(X, gamma, kernel = "h_w") {
    Xv <- if (kernel == "h_0") abs(X) else X
    n  <- length(Xv)
    H  <- outer(Xv, Xv, function(a, b) sign(a - b))
    sigma <- sqrt(mean(colMeans(H)^2))
    w  <- holder_weights(n, gamma, 0, symmetric = TRUE)
    best <- 0
    for (k in 0:(n - 1L)) for (m in (k + 1L):n) {
      I <- (k + 1L):m
      v <- abs(sum(H[I, setdiff(1:n, I)])) * w[m - k]
      if (is.finite(v) && v > best) best <- v
    }
    best / (n^1.5 * sigma)
  }
  set.seed(SEED); x <- rnorm(60)
  for (g in c(0, 0.25, 0.45)) {
    expect_equal(wilcoxon_statistic(x, g), brute_wilcoxon(x, g),
                 tolerance = TOL)
  }
  expect_equal(wilcoxon_statistic(x, 0.25, kernel = "h_0"),
               brute_wilcoxon(x, 0.25, kernel = "h_0"), tolerance = TOL)
  # vectorised gamma and the scan agree with the scalar statistic
  gb <- wilcoxon_statistic(x, c(0, 0.25, 0.45))
  st <- vapply(c(0, 0.25, 0.45), function(g) wilcoxon_statistic(x, g), 0)
  sc <- vapply(c(0, 0.25, 0.45), function(g) wilcoxon_scan(x, g)$stat, 0)
  expect_equal(gb, st, tolerance = TOL)
  expect_equal(sc, st, tolerance = TOL)
})

test_that("lp_statistic equals the printed centred formula, brute-forced", {
  brute_lp <- function(X, gamma, p = 2, cdf = pnorm, t_grid = 100) {
    n  <- length(X); F.X <- cdf(X)
    w  <- holder_weights(n, gamma)
    ts <- seq(0, 1, length.out = t_grid + 1)
    ind <- outer(F.X, ts, function(f, t) (f <= t) - t)
    tot <- colSums(ind)
    best <- 0
    for (k in 0:(n - 1L)) for (l in 1:(n - k)) {
      cs <- colSums(ind[(k + 1L):(k + l), , drop = FALSE])
      cs <- cs - (l / n) * tot
      v <- (sum(abs(cs)^p) / t_grid)^(1 / p) * w[l]
      if (is.finite(v) && v > best) best <- v
    }
    best / sqrt(n)
  }
  set.seed(SEED); x <- rnorm(60)
  for (g in c(0, 0.25, 0.45)) {
    expect_equal(lp_statistic(x, g, cdf = pnorm), brute_lp(x, g),
                 tolerance = TOL)
  }
  expect_equal(lp_statistic(x, 0.25, p = 4, cdf = pnorm),
               brute_lp(x, 0.25, p = 4), tolerance = TOL)
  expect_equal(lp_scan(x, 0.25, cdf = pnorm)$stat,
               lp_statistic(x, 0.25, cdf = pnorm), tolerance = TOL)
})

test_that("rkhs_statistic equals the printed formula, brute-forced", {
  brute_rkhs <- function(X, K, gamma) {
    n <- length(X)
    w <- holder_weights(n, gamma)
    best <- 0
    for (l in 1:(n %/% 2L)) for (a in 1:(n - l + 1L)) {
      I  <- a:(a + l - 1L)
      ln <- l / n
      d2 <- sum(K[I, I]) + ln^2 * sum(K) - 2 * ln * sum(K[I, 1:n])
      v  <- sqrt(max(d2, 0)) * w[l]
      if (is.finite(v) && v > best) best <- v
    }
    best / sqrt(n)
  }
  set.seed(SEED); x <- rnorm(60)
  Kg <- gram(x, kernel_gaussian)
  for (g in c(0, 0.25, 0.45)) {
    expect_equal(rkhs_statistic(x, g, kernel = "gaussian"),
                 brute_rkhs(x, Kg, g), tolerance = TOL)
  }
  K2 <- gram(x, kernel_2max)
  expect_equal(rkhs_statistic(x, 0.25, kernel = "2max"),
               brute_rkhs(x, K2, 0.25), tolerance = TOL)
  expect_equal(rkhs_scan(x, 0.25)$stat,
               rkhs_statistic(x, 0.25), tolerance = TOL)
})

test_that("far1_statistic equals the printed formula, brute-forced", {
  brute_far1 <- function(X, gamma, beta_t = far1_betas$beta1) {
    L <- nrow(X); n <- ncol(X)
    s <- seq(0, 1, length.out = L)
    W <- (1 - beta_t(s)) * (X - rowMeans(X))
    best <- 0
    for (k in 0:(n - 1L)) for (m in (k + 1L):n) {
      Ssum  <- rowSums(W[, (k + 1L):m, drop = FALSE])
      v <- max(abs(Ssum)) / (m - k)^gamma
      if (v > best) best <- v
    }
    best * n^(-0.5 + gamma)
  }
  set.seed(SEED); Xf <- generate_f4(40, 0, c(0.4, 0.6))
  for (g in c(0.25, 0.4)) {
    expect_equal(far1_statistic(Xf, g), brute_far1(Xf, g),
                 tolerance = TOL)
  }
  expect_equal(far1_scan(Xf, 0.4)$stat, far1_statistic(Xf, 0.4),
               tolerance = TOL)
})
