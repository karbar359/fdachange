test_that("method = none returns the input unchanged", {
  set.seed(1)
  X <- generate_f1(10, a = 0, l = 20)
  R <- reconstruct(X, method = "none")
  expect_identical(R$X, X)
  expect_equal(R$argvals, seq(0, 1, length.out = 21))
})

test_that("fs_interp is a fixed point on polygonal input", {
  # A curve already piecewise-linear on the dyadic mesh 2^-K is reproduced
  # exactly by the level-K Faber-Schauder truncation.
  K <- 4L
  grid <- seq(0, 1, by = 2^(-K))          # 17 vertices
  set.seed(2)
  V <- matrix(rnorm(length(grid) * 8), nrow = length(grid))
  R <- reconstruct(V, method = "fs_interp", K = K, grid = "vertices",
                   argvals = grid)
  expect_equal(R$X, V, tolerance = 1e-12)
  # and interpolation back onto a finer input grid stays polygonal
  fine <- seq(0, 1, length.out = 65)
  Vf <- apply(V, 2L, function(col) approx(grid, col, xout = fine)$y)
  Rf <- reconstruct(Vf, method = "fs_interp", K = K, grid = "vertices",
                    argvals = fine)
  expect_equal(Rf$X, V, tolerance = 1e-12)
})

test_that("fs_interp reduces to the dyadic vertex grid", {
  set.seed(3)
  X <- generate_f1(6, a = 0, l = 100)
  R <- reconstruct(X, method = "fs_interp", K = 4)
  expect_equal(dim(R$X), c(17L, 6L))
  expect_equal(R$argvals, seq(0, 1, by = 2^(-4)))
  Ri <- reconstruct(X, method = "fs_interp", K = 4, grid = "input")
  expect_equal(dim(Ri$X), dim(X))
})

test_that("fda-backed reconstructions run when fda is present", {
  skip_if_not_installed("fda")
  set.seed(4)
  X <- generate_f1(6, a = 0, l = 100)
  R <- reconstruct(X, method = "fs_ls", K = 4)
  expect_equal(dim(R$X), c(17L, 6L))
  Rp <- reconstruct(X, method = "bspline_pen", K = 4, lambda = 1e-4)
  expect_equal(dim(Rp$X), c(17L, 6L))
  # least squares reproduces polygonal input on its own mesh
  grid <- seq(0, 1, by = 2^(-4))
  V <- matrix(rnorm(17 * 5), nrow = 17)
  Rv <- reconstruct(V, method = "fs_ls", K = 4, argvals = grid)
  expect_equal(Rv$X, V, tolerance = 1e-6)
})
