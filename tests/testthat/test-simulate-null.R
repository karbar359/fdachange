test_that("simulators return the requested number of draws", {
  set.seed(1)
  d <- simulate_wilcoxon_null(5, gamma = 0.25, grid = 100, progress = FALSE)
  expect_length(d, 5)
  expect_true(all(is.finite(d)) && all(d > 0))
  d <- simulate_lp_null(3, gamma = 0.25, p = 2, grid = 30, progress = FALSE)
  expect_length(d, 3)
  d <- simulate_rkhs_null(3, gamma = 0.25, kernel = "gaussian", n = 80,
                          progress = FALSE)
  expect_length(d, 3)
  d <- simulate_birbilas_null(3, gamma = 0.4, grid = 40, progress = FALSE)
  expect_length(d, 3)
})

test_that("simulators are reproducible under set.seed", {
  set.seed(42); a <- simulate_wilcoxon_null(3, 0.25, grid = 50,
                                            progress = FALSE)
  set.seed(42); b <- simulate_wilcoxon_null(3, 0.25, grid = 50,
                                            progress = FALSE)
  expect_identical(a, b)
})

test_that("wilcoxon simulator agrees with the shipped cache distribution", {
  skip_if(!nzchar(Sys.getenv("FDACHANGE_RUN_SLOW")),
          "slow: set FDACHANGE_RUN_SLOW=1")
  set.seed(7)
  d <- simulate_wilcoxon_null(400, gamma = 0.25, grid = 2000,
                              progress = FALSE)
  ks <- suppressWarnings(stats::ks.test(d, wilcoxon_null_draws(0.25)))
  expect_gt(ks$p.value, 0.001)
})
