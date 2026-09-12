# Curve reconstruction (smoothing) as an explicit, parameterised pipeline
# step. The level-K Faber-Schauder truncation and the span of order-2
# B-splines on the dyadic mesh 2^-K name the SAME function space
# (continuous piecewise-linear functions, dimension 2^K + 1) but different
# maps into it: FS truncation interpolates at the dyadic points
# (fs_interp), Data2fd least-squares over the whole sampling grid (fs_ls).
# On rough curves the two differ materially.

.recon_methods <- c("none", "fs_interp", "fs_ls", "bspline_pen")
.recon_lambda_grid <- 10^seq(-8, 2, by = 1)

.recon_env <- new.env(parent = emptyenv())

# Order-2 (piecewise-linear) B-spline basis on the dyadic mesh 2^-K;
# dimension 2^K + 1. Cached per K.
.fs_basis <- function(K) {
  key <- paste0("fs_K", K)
  if (is.null(.recon_env[[key]])) {
    .recon_env[[key]] <- fda::create.bspline.basis(
      norder = 2L, breaks = seq(0, 1, by = 2^(-as.integer(K))))
  }
  .recon_env[[key]]
}

# Order-4 (cubic) B-spline basis on the same mesh, for the penalised arm:
# the only genuinely smooth reconstruction (order 2 has no curvature to
# penalise).
.pen_basis <- function(K) {
  key <- paste0("pen_K", K)
  if (is.null(.recon_env[[key]])) {
    .recon_env[[key]] <- fda::create.bspline.basis(
      norder = 4L, breaks = seq(0, 1, by = 2^(-as.integer(K))))
  }
  .recon_env[[key]]
}

.recon_grid <- function(method, K, grid, n_input) {
  if (identical(method, "none") || identical(grid, "input")) {
    seq(0, 1, length.out = n_input)
  } else {
    seq(0, 1, by = 2^(-as.integer(K)))
  }
}

# Pick lambda by pooled GCV over .recon_lambda_grid.
.gcv_lambda <- function(argvals, X, basis) {
  scores <- vapply(.recon_lambda_grid, function(lam) {
    fdp <- fda::fdPar(basis, fda::int2Lfd(2L), lam)
    sum(fda::smooth.basis(argvals, X, fdp)$gcv, na.rm = TRUE)
  }, numeric(1))
  .recon_lambda_grid[which.min(scores)]
}

#' Reconstruct (smooth) a sample of curves
#'
#' Applies one deterministic, curve-wise reconstruction to every
#' observation before projection. Reconstruction preserves the i.i.d.
#' structure and the location of a changed segment; it changes which
#' distribution is being tested, not whether the epidemic framework
#' applies.
#'
#' Methods:
#' * `"none"`: return the sampled curves unchanged.
#' * `"fs_interp"`: level-`K` Faber-Schauder truncation, i.e. polygonal
#'   interpolation at the dyadic points `k / 2^K`. Base R.
#' * `"fs_ls"`: least squares onto the same level-`K` polygonal space
#'   (order-2 B-splines via `fda::Data2fd`); this is what produced the
#'   dissertation's application tables. Needs fda.
#' * `"bspline_pen"`: order-4 B-splines with a second-derivative roughness
#'   penalty; the only genuinely smooth arm. Needs fda.
#'
#' The evaluation grid matters: `grid = "vertices"` evaluates on the
#' `2^K + 1` dyadic points (reproduces the printed application tables),
#' `grid = "input"` evaluates back on the original sampling grid.
#'
#' @param X numeric matrix, one curve per **column**, rows = sampling
#'   grid.
#' @param method one of `"none"`, `"fs_interp"`, `"fs_ls"`,
#'   `"bspline_pen"`.
#' @param K dyadic truncation level; the reconstruction space has
#'   dimension `2^K + 1`.
#' @param grid `"vertices"` or `"input"`, see Details.
#' @param lambda roughness penalty for `"bspline_pen"`; `NULL` selects by
#'   GCV over `10^(-8..2)`.
#' @param argvals sampling grid of the input rows (default: equispaced on
#'   `[0, 1]`).
#' @return a list with `X` (matrix on the evaluation grid), `argvals`
#'   (that grid), `method`, `K`, `grid`.
#' @references Bartkus (2026), Chapters 1 and 7.
#' @examples
#' set.seed(1)
#' X <- generate_f1(20, a = 0)
#' R <- reconstruct(X, method = "fs_interp", K = 4)
#' dim(R$X)   # 17 x 20: curves on the dyadic vertex grid
#' @export
reconstruct <- function(X, method = c("none", "fs_interp", "fs_ls",
                                      "bspline_pen"),
                        K = 4L, grid = c("vertices", "input"), lambda = 1e-4,
                        argvals = NULL) {
  method <- match.arg(method)
  grid <- match.arg(grid)
  X <- as.matrix(X)
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = nrow(X))
  out_grid <- .recon_grid(method, K, grid, nrow(X))

  if (identical(method, "none")) {
    return(list(X = X, argvals = argvals, method = method, K = NA_integer_,
                grid = "input"))
  }
  if (method %in% c("fs_ls", "bspline_pen")) {
    if (!requireNamespace("fda", quietly = TRUE)) {
      stop("reconstruct(method = \"", method, "\") needs the 'fda' package",
           call. = FALSE)
    }
  }

  Xr <- switch(method,
    fs_interp = {
      # FS truncation == polygonal interpolation at the dyadic vertices.
      V <- apply(X, 2L, function(col)
        approx(argvals, col, xout = seq(0, 1, by = 2^(-as.integer(K))))$y)
      V <- matrix(V, nrow = 2^as.integer(K) + 1L)
      if (identical(grid, "input")) {
        apply(V, 2L, function(col)
          approx(seq(0, 1, by = 2^(-as.integer(K))), col, xout = argvals)$y)
      } else {
        V
      }
    },
    fs_ls = {
      fd <- fda::Data2fd(argvals, X, basisobj = .fs_basis(K))
      fda::eval.fd(out_grid, fd)
    },
    bspline_pen = {
      b <- .pen_basis(K)
      lam <- if (is.null(lambda)) .gcv_lambda(argvals, X, b) else lambda
      fdp <- fda::fdPar(b, fda::int2Lfd(2L), lam)
      fd  <- fda::smooth.basis(argvals, X, fdp)$fd
      fda::eval.fd(out_grid, fd)
    })

  Xr <- matrix(as.numeric(Xr), nrow = length(out_grid))
  list(X = Xr, argvals = out_grid, method = method, K = as.integer(K),
       grid = grid)
}
