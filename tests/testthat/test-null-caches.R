test_that("shipped caches load and expose the documented cells", {
  expect_length(wilcoxon_null_draws(0.25), 30000L)
  expect_length(wilcoxon_null_draws(0.5, beta = 1), 30000L)
  expect_length(lp_null_draws(0.45, p = 10), 300L)
  expect_length(rkhs_null_draws(0.45, "2max", "rbeta.1.9"), 1000L)
  expect_length(rkhs_null_draws(0, "min", "runif"), 1000L)
  expect_length(birbilas_null_draws(0.4), 1000L)
})

test_that("gamma = 0 Wilcoxon null is the Kuiper grid, not the cache", {
  d <- wilcoxon_null_draws(0)
  expect_length(d, 200000L)
  expect_false(is.unsorted(d))   # deterministic quantile grid
})

test_that("a missing cell errors and names the simulator", {
  expect_error(wilcoxon_null_draws(0.33), "simulate_wilcoxon_null")
  expect_error(lp_null_draws(0.25, p = 3), "simulate_lp_null")
  expect_error(rkhs_null_draws(0.25, "min"), "simulate_rkhs_null")
  expect_error(rkhs_null_draws(0.25, "gaussian", "runif"),
               "simulate_rkhs_null")
  expect_error(birbilas_null_draws(0.3), "simulate_birbilas_null")
})

test_that("save_null_draws stores into the user cache and takes precedence", {
  cache <- file.path(tempdir(), "fdachange-test-cache")
  old_opt <- options(fdachange.cache_dir = cache)
  on.exit({
    options(old_opt)
    unlink(cache, recursive = TRUE)
    rm(list = ls(envir = fdachange:::.null_env), envir = fdachange:::.null_env)
  })

  expect_equal(null_cache_dir(), cache)
  # new cell becomes available
  set.seed(1); draws <- runif(500)
  save_null_draws(draws, "wilcoxon", gamma = 0.33)
  expect_equal(wilcoxon_null_draws(0.33), draws)
  # user entry overrides a shipped cell
  save_null_draws(1:100 / 100, "lp", gamma = 0.25, p = 2)
  expect_length(lp_null_draws(0.25, p = 2), 100L)
  # rkhs keying round-trips through kernel + baseline
  save_null_draws(draws, "rkhs", gamma = 0.1, kernel = "2max",
                  baseline = "runif")
  expect_equal(rkhs_null_draws(0.1, "2max", "runif"), draws)
})

test_that("null draws feed critical values consistently", {
  d <- lp_null_draws(0.25, p = 2)
  expect_equal(lp_critical_value(0.25, p = 2),
               as.numeric(quantile(d, 0.95)))
  expect_equal(lp_critical_value(0.25, p = 2, splits = 2),
               as.numeric(quantile(d, sqrt(0.95))))
})
