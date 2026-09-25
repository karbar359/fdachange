test_that("epidemic_test returns a coherent htest object", {
  set.seed(3)
  x <- rnorm(200); x[81:140] <- x[81:140] + 1
  for (m in c("wilcoxon", "lp", "rkhs")) {
    ht <- epidemic_test(x, m, gamma = 0.25)
    expect_s3_class(ht, c("epidemic_htest", "htest"))
    expect_named(ht$statistic, "T")
    expect_true(ht$p.value >= 0 && ht$p.value <= 1)
    # p < alpha exactly when the statistic clears the critical value
    expect_equal(ht$p.value < ht$alpha,
                 unname(ht$statistic) > ht$critical.value)
    # the detected window is inside the sample
    expect_true(ht$k >= 0 && ht$m <= 200 && ht$k < ht$m)
  }
})

test_that("a strong epidemic is detected and localised", {
  set.seed(4)
  x <- rnorm(200); x[81:140] <- x[81:140] + 3
  ht <- epidemic_test(x, "wilcoxon", gamma = 0.25)
  expect_lt(ht$p.value, 0.01)
  expect_gt(ht$m, 100)          # window overlaps the true segment
  expect_lt(ht$k, 120)
})

test_that("wilcoxon gamma = 0 uses the exact Kuiper law", {
  set.seed(3)
  x <- rnorm(150)
  ht <- epidemic_test(x, "wilcoxon", gamma = 0)
  expect_identical(ht$n.null, Inf)
  expect_equal(1 - ht$p.value, kuiper_cdf(unname(ht$statistic)),
               tolerance = 1e-12)
  expect_equal(ht$critical.value, kuiper_quantile(0.95), tolerance = 1e-9)
})

test_that("null = injection overrides the automatic calibration", {
  set.seed(5)
  x <- rnorm(100)
  draws <- seq(0, 10, length.out = 1001)
  ht <- epidemic_test(x, "rkhs", gamma = 0.25, null = draws)
  expect_equal(ht$critical.value, as.numeric(quantile(draws, 0.95)))
  expect_equal(ht$n.null, 1001L)
  expect_equal(1 - ht$p.value, mean(draws <= unname(ht$statistic)))
})

test_that("splits > 1 aggregates chunk maxima and adjusts the calibration", {
  set.seed(6)
  x <- rnorm(200)
  ht2 <- epidemic_test(x, "wilcoxon", gamma = 0.25, splits = 2)
  s1 <- wilcoxon_statistic(x[1:100], 0.25)
  s2 <- wilcoxon_statistic(x[101:200], 0.25)
  expect_equal(unname(ht2$statistic), max(s1, s2))
  expect_equal(ht2$critical.value,
               wilcoxon_critical_value(0.25, splits = 2))
  expect_true(is.na(ht2$k))     # no single window is reported
  d <- wilcoxon_null_draws(0.25)
  expect_equal(ht2$p.value, 1 - mean(d <= max(s1, s2))^2)
})

test_that("far1 method consumes the functional matrix directly", {
  set.seed(7)
  X <- generate_f4(60, a = 2, l = 50)
  ht <- epidemic_test(X, "far1", gamma = 0.4)
  expect_equal(unname(ht$statistic), far1_statistic(X, 0.4))
  expect_true(ht$k >= 0 && ht$m <= 60)
})

test_that("rkhs with a custom kernel function requires explicit null draws", {
  set.seed(8)
  x <- rnorm(50)
  my_kernel <- function(a, b) exp(-abs(a - b))
  expect_error(epidemic_test(x, "rkhs", gamma = 0.25, kernel = my_kernel),
               "null")
  ht <- epidemic_test(x, "rkhs", gamma = 0.25, kernel = my_kernel,
                      null = runif(200, 0, 5))
  expect_s3_class(ht, "htest")
})

test_that("print method shows the critical value and the window", {
  set.seed(9)
  x <- rnorm(100); x[41:70] <- x[41:70] + 2
  ht <- epidemic_test(x, "wilcoxon", gamma = 0.25)
  out <- capture.output(print(ht))
  expect_true(any(grepl("critical value", out)))
  expect_true(any(grepl("detected window", out)))
})
