test_that("Kuiper closed form matches the dissertation's Table 3.1 row", {
  printed_levels <- c(0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999)
  printed_row    <- c(1.213, 1.460, 1.612, 1.741, 1.857, 2.012, 2.104, 2.195,
                      2.311)
  qk <- kuiper_quantile(printed_levels)
  # Bulk levels: grid discretisation of the sup biases the printed MC row
  # down, so printed <= exact there; the extreme tail is MC noise. Both
  # regimes stay within 0.014.
  expect_true(all(abs(printed_row - qk) < 0.014))
  expect_true(all((printed_row - qk)[1:5] <= 0))
})

test_that("the quantile-grid stand-in reproduces the exact quantiles", {
  g <- wilcoxon_null_draws(gamma = 0)
  expect_equal(as.numeric(quantile(g, 0.95)), 1.7473, tolerance = 5e-4)
  expect_equal(as.numeric(quantile(g, 0.5)), kuiper_quantile(0.5),
               tolerance = 1e-3)
})

test_that("cdf and quantile invert each other", {
  v <- c(1.2, 1.7473, 2.0)
  expect_equal(kuiper_quantile(kuiper_cdf(v)), v, tolerance = 1e-6)
  expect_equal(kuiper_cdf(0.1), 0)   # below support
})
