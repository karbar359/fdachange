test_that("scalar generators return length n; functional (l+1) x n", {
  set.seed(1)
  expect_length(generate_s1(50), 50)
  expect_length(generate_s2(50, a = 1), 50)
  expect_length(generate_s3(50, a = 1), 50)
  expect_length(generate_s4(50, a = 1), 50)
  expect_length(generate_s5(50, a = 1), 50)
  expect_equal(dim(generate_f1(20, a = 3, l = 30)), c(31L, 20L))
  expect_equal(dim(generate_f2(20, a = 0.5, l = 30)), c(31L, 20L))
  expect_equal(dim(generate_f3(20, a = 3, l = 30)), c(31L, 20L))
  expect_equal(dim(generate_f4(20, a = 2, l = 30)), c(31L, 20L))
  expect_equal(dim(generate_f5(20, a = 2, l = 30)), c(31L, 20L))
})

# F5 has no golden-value entries: the archive's Wiener + spike scenario was
# never part of the audited pipeline, so these tests pin the restored
# implementation directly.
test_that("F5 adds the spike to window curves only, at mid-curve height a", {
  set.seed(5)
  n <- 400; a <- 5
  X <- generate_f5(n, a = a, ec_interval = c(0.4, 0.7))
  mid <- 51L                                  # t = 0.5 on the (0:100)/100 grid
  win <- 161:280                              # exact window {k*+1..k*+l}
  out <- setdiff(seq_len(n), win)
  # spike(0.5) = 1, Wiener sd at t = 0.5 is sqrt(0.5); averaging over the
  # window columns leaves the +a shift visible well beyond MC noise
  expect_equal(mean(X[mid, win]) - mean(X[mid, out]), a, tolerance = 0.1)
  # non-window curves start at 0 (standard Wiener); window curves start at
  # a * f_alpha(0) because the dissertation's spike has nonzero shoulders
  expect_true(all(X[1, out] == 0))
  f0 <- (1 - 0.5^0.45)^(1 / 0.45)
  expect_equal(unname(X[1, win]), rep(a * f0, length(win)))
})

test_that("F5 spike alpha changes only the deterministic bump (same stream)", {
  s <- (0:100) / 100
  f_alpha <- function(t, al) (1 - abs(t - 0.5)^al)^(1 / al)
  set.seed(9); Xa <- generate_f5(50, a = 3, spike_alpha = 0.45)
  set.seed(9); Xb <- generate_f5(50, a = 3, spike_alpha = 0.2)
  diff_expected <- 3 * (f_alpha(s, 0.45) - f_alpha(s, 0.2))
  win <- 21:35                                # exact window at n = 50
  for (j in win) expect_equal(Xa[, j] - Xb[, j], diff_expected)
  out <- setdiff(1:50, win)
  expect_true(all(Xa[, out] == Xb[, out]))
  # the spike is the dissertation's f_alpha: unit peak at t = 1/2
  expect_equal(f_alpha(0.5, 0.45), 1)
})

test_that("F5 exact window holds l columns, legacy l + 1", {
  changed_cols <- function(exact) {
    set.seed(1)
    X <- generate_f5(400L, a = 50, ec_interval = c(0.4925, 0.5075),
                     exact_window = exact)
    sum(X[51L, ] > 25)                        # spike height 50 at t = 0.5
  }
  expect_equal(changed_cols(TRUE), 6L)
  expect_equal(changed_cols(FALSE), 7L)
})

test_that("exact window holds l points, legacy window l + 1", {
  # Short regime: l = 6 at n = 400 with the interval (0.4925, 0.5075).
  count_changed <- function(exact) {
    set.seed(1)
    x <- generate_s1(400L, 50, c(0.4925, 0.5075), exact_window = exact)
    sum(x > 25)
  }
  expect_equal(count_changed(FALSE), 7L)
  expect_equal(count_changed(TRUE), 6L)
})

test_that("S2 generates noncentrality a^2 on the underlying normal", {
  set.seed(7)
  a <- 2
  x <- generate_s2(200000L, a, c(0.1, 0.9))
  m <- mean(x[20001:180000])
  expect_equal(m, 1 + a^2, tolerance = 0.05)
})

test_that("scenario_cdf returns the matching null CDF", {
  expect_equal(scenario_cdf("S1")(c(-1, 0, 1)), pnorm(c(-1, 0, 1)))
  expect_equal(scenario_cdf("S2")(c(0.5, 2)), pchisq(c(0.5, 2), df = 1))
  # Pareto(1, 3): F(x) = 1 - x^{-3} for x >= 1
  expect_equal(scenario_cdf("S3")(c(1.5, 3)), 1 - c(1.5, 3)^(-3))
  expect_equal(scenario_cdf("S3")(0.5), 0)
  expect_equal(scenario_cdf("S5")(2), scenario_cdf("S3")(2))
  expect_equal(scenario_cdf("S4")(0.2), pbeta(0.2, 1, 9))
  x <- rnorm(20)
  expect_equal(scenario_cdf("F1")(x), ecdf(x)(x))
  expect_equal(scenario_cdf("F5")(x), ecdf(x)(x))
  expect_error(scenario_cdf("X9"), "unknown scenario")
})

test_that("the change window actually changes the distribution", {
  set.seed(2)
  x <- generate_s1(300, a = 5, ec_interval = c(0.3, 0.6))
  expect_gt(mean(x[92:180]), mean(x[c(1:90, 182:300)]) + 3)
})
