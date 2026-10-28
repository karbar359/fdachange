# The parallel helper gives every item its own L'Ecuyer substream, so draws
# are independent, identical for any core count, and reproducible under
# set.seed() whatever generator the caller uses.

draw_one <- function(items, cores) {
  unlist(fdachange:::.apply_progress(items, function(i) rnorm(1),
                                     cores = cores))
}

test_that("every item draws from its own substream", {
  skip_on_cran()
  skip_on_os("windows")
  set.seed(2)
  x <- draw_one(1:24, cores = 2)
  expect_length(unique(x), 24L)
})

test_that("results do not depend on the core count", {
  skip_on_cran()
  skip_on_os("windows")
  set.seed(1); a <- draw_one(1:24, cores = 1)
  set.seed(1); b <- draw_one(1:24, cores = 2)
  expect_identical(a, b)
})

test_that("the caller's generator kind and stream survive the call", {
  old <- RNGkind()
  on.exit(RNGkind(old[1L], old[2L], old[3L]))
  RNGkind("Mersenne-Twister")
  set.seed(5); a <- draw_one(1:6, cores = 1); after_a <- runif(1)
  set.seed(5); b <- draw_one(1:6, cores = 1); after_b <- runif(1)
  expect_identical(RNGkind()[1L], "Mersenne-Twister")
  expect_identical(a, b)
  expect_identical(after_a, after_b)
})

test_that("an L'Ecuyer caller is left on a fresh stream", {
  old <- RNGkind()
  on.exit(RNGkind(old[1L], old[2L], old[3L]))
  RNGkind("L'Ecuyer-CMRG"); set.seed(3)
  a <- draw_one(1:4, cores = 1)
  b <- draw_one(1:4, cores = 1)
  expect_false(any(a %in% b))
})
