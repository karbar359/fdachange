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
  expect_error(scenario_cdf("X9"), "unknown scenario")
})

test_that("the change window actually changes the distribution", {
  set.seed(2)
  x <- generate_s1(300, a = 5, ec_interval = c(0.3, 0.6))
  expect_gt(mean(x[92:180]), mean(x[c(1:90, 182:300)]) + 3)
})
