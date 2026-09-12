test_that("base-R projections have the right shape and determinism", {
  set.seed(1)
  X <- generate_f1(30, a = 0, l = 20)
  l2 <- project_l2(X)
  expect_length(l2, 30)
  expect_true(all(l2 >= 0))
  expect_identical(project_l2(X), l2)     # deterministic
  set.seed(2); i1 <- project_ito(X)
  set.seed(2); i2 <- project_ito(X)
  expect_identical(i1, i2)                # same path, same projection
  set.seed(3)
  expect_false(identical(project_ito(X), i1))  # fresh random direction
})

test_that("projection_fn resolves ids and multiplicity is right", {
  expect_identical(projection_fn("l2"), project_l2)
  expect_identical(projection_fn("none")(1:5), 1:5)
  expect_error(projection_fn("nope"), "unknown projection")
  expect_equal(projection_multiplicity("l2"), 1L)
  expect_equal(projection_multiplicity("pc_d", 4L), 4L)
})

test_that("fda-backed projections work when fda.usc is present", {
  skip_if_not_installed("fda")
  skip_if_not_installed("fda.usc")
  set.seed(4)
  X <- generate_f1(25, a = 0, l = 20)
  expect_length(project_c_max_var(X), 25)
  expect_length(project_wiener(X), 25)
  expect_length(project_pc1(X), 25)
  pd <- project_pc_d(X, d = 2)
  expect_equal(dim(pd), c(2L, 25L))
})

test_that("make_run and scan_over_splits compose projections and statistics", {
  set.seed(5)
  X <- generate_f1(60, a = 0, l = 20)
  run <- make_run(project_l2, rkhs_statistic, gamma = 0.25,
                  stat_args = list(kernel = "gaussian"))
  expect_equal(scan_over_splits(X, 1, run), run(X))
  v2 <- scan_over_splits(X, 2, run)
  expect_equal(v2, max(run(X[, 1:30]), run(X[, 31:60])))
  expect_error(scan_over_splits(X, 7, run), "not divisible")
  # matrix projection: statistic per row, max over rows
  fake_proj <- function(Z) rbind(project_l2(Z), 2 * project_l2(Z))
  run_m <- make_run(fake_proj, rkhs_statistic, gamma = 0.25,
                    stat_args = list(kernel = "gaussian"))
  expect_equal(run_m(X),
               max(rkhs_statistic(project_l2(X), 0.25),
                   rkhs_statistic(2 * project_l2(X), 0.25)))
  # matrix_in forwards the whole matrix
  run_b <- make_run(identity, birbilas_statistic, gamma = 0.4,
                    matrix_in = TRUE)
  expect_equal(run_b(X), birbilas_statistic(X, 0.4))
})
