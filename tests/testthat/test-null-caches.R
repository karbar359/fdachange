test_that("shipped caches load and expose the documented cells", {
  expect_length(wilcoxon_null_draws(0.25), 30000L)
  expect_length(wilcoxon_null_draws(0.5, beta = 1), 30000L)
  expect_length(lp_null_draws(0.45, p = 10), 300L)
  expect_length(rkhs_null_draws(0.45, "2max", "rbeta.1.9"), 1000L)
  expect_length(rkhs_null_draws(0, "min", "runif"), 1000L)
  expect_length(far1_null_draws(0.4), 1000L)
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
  expect_error(far1_null_draws(0.3), "simulate_far1_null")
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

test_that("shipped cache cells are independent draws", {
  # Reused RNG substreams show up as autocorrelation at the batch lag: the
  # pre-fix RKHS and far1 caches had 0.84 to 0.92 (about 100 distinct draws
  # of 1000). The gamma = 0 Wilcoxon cell is the Kuiper grid, not a cache.
  for (m in c("wilcoxon", "lp", "rkhs", "far1")) {
    tbl <- readRDS(system.file("extdata", paste0(m, "_null.rds"),
                               package = "fdachange"))
    for (k in names(tbl)) {
      x <- tbl[[k]]
      a <- stats::acf(x, lag.max = min(200L, length(x) %/% 3L),
                      plot = FALSE)$acf[-1]
      expect_lt(max(abs(a)), 0.25, label = sprintf("%s %s max |acf|", m, k))
    }
  }
  wilcoxon <- readRDS(system.file("extdata", "wilcoxon_null.rds",
                                  package = "fdachange"))
  expect_false("g0_b0" %in% names(wilcoxon))
})
