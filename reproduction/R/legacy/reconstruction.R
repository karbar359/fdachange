#' Curve reconstruction (smoothing) as an explicit, parameterised pipeline step.
#'
#' The dissertation describes the preprocessing step as a Faber-Schauder basis
#' expansion truncated at level K (Ch. 1 "Smoothing with Faber-Schauder basis";
#' Ch. 7 Sec. 7.2, K = 4). The original code implements it as
#' `create.bspline.basis(norder = 2, breaks = seq(0, 1, by = 2^-K))` fed to
#' `Data2fd`. Those two descriptions name the SAME function space:
#'
#'   span{ order-2 B-splines, breaks k/2^K }  ==  span{ FS system truncated at K }
#'                                            ==  continuous piecewise-linear
#'                                                functions on the mesh 2^-K,
#'   both of dimension 2^K + 1.
#'
#' They do NOT name the same map into that space. Faber-Schauder truncation is
#' interpolation at the dyadic points; `Data2fd` is least squares over the whole
#' sampling grid. On rough curves the two differ materially. Hence the split
#' between `fs_interp` and `fs_ls` below.
#'
#' A third axis matters and is easy to miss: the EVALUATION grid. A polygonal
#' function is determined by its vertices, but the trapezoidal L2 norm over
#' 2^K + 1 vertices is not the same number as the trapezoidal norm over the
#' original sampling grid (about 7% apart on the Litgrid weekly curves). The
#' printed application tables use the vertex grid, which is what
#' `fda.usc::fdata()` produces by default from a 2^K + 1 dimensional basis.
#'
#' Reconstruction is a deterministic, curve-wise map applied identically to
#' every observation, so it preserves the i.i.d. structure and the location of
#' any changed segment. It changes which distribution is being tested, not
#' whether the epidemic hypothesis framework applies.
#'
#' Callers reconstruct FIRST and then project; `shared/projections.R` consumes
#' whatever grid it is handed.

suppressPackageStartupMessages({
  library(fda)
})

# ---- Defaults (edit these) ------------------------------------------------
RECON_METHOD_DEFAULT <- "none"   # none | fs_interp | fs_ls | bspline_pen
RECON_K_DEFAULT      <- 4L       # dyadic truncation level; space dim = 2^K + 1
RECON_GRID_DEFAULT   <- "vertices"  # vertices | input
RECON_LAMBDA_DEFAULT <- 1e-4     # bspline_pen only; NULL selects by GCV
# ---------------------------------------------------------------------------

RECONSTRUCTION_METHODS <- c("none", "fs_interp", "fs_ls", "bspline_pen")

#' Order-2 (piecewise-linear) B-spline basis on the dyadic mesh 2^-K.
#' Cached per K. `nbasis` is 2^K + 1, the dimension of the level-K
#' Faber-Schauder truncation.
.fs_basis_cache <- new.env(parent = emptyenv())
fs_basis <- function(K) {
  key <- paste0("K", K)
  if (is.null(.fs_basis_cache[[key]])) {
    .fs_basis_cache[[key]] <- create.bspline.basis(
      norder = 2L, breaks = seq(0, 1, by = 2^(-as.integer(K))))
  }
  .fs_basis_cache[[key]]
}

#' Order-4 (cubic) B-spline basis on the same dyadic mesh, for the penalised
#' arm. This is the only genuinely smooth reconstruction in the pipeline: with
#' order 2 the second derivative vanishes identically, so a curvature penalty
#' is vacuous no matter how lambda is set.
.pen_basis_cache <- new.env(parent = emptyenv())
pen_basis <- function(K) {
  key <- paste0("K", K)
  if (is.null(.pen_basis_cache[[key]])) {
    .pen_basis_cache[[key]] <- create.bspline.basis(
      norder = 4L, breaks = seq(0, 1, by = 2^(-as.integer(K))))
  }
  .pen_basis_cache[[key]]
}

#' The evaluation grid for a given method/level/grid choice.
recon_grid <- function(method, K, grid, n_input) {
  if (identical(method, "none") || identical(grid, "input")) {
    seq(0, 1, length.out = n_input)
  } else {
    seq(0, 1, by = 2^(-as.integer(K)))
  }
}

#' Reconstruct a sample of curves.
#'
#' @param X matrix, one curve per COLUMN, rows = sampling grid (the pipeline's
#'   functional-scenario contract).
#' @param method one of `RECONSTRUCTION_METHODS`:
#'   - `none`        raw sampled curves, returned unchanged
#'   - `fs_interp`   Faber-Schauder truncation: piecewise-linear interpolation
#'                   at the dyadic points k/2^K. This is what the thesis text
#'                   describes.
#'   - `fs_ls`       least squares onto the same level-K polygonal space via
#'                   `Data2fd`. This is what the original code does, and what
#'                   produced the printed application tables.
#'   - `bspline_pen` order-4 B-splines with a second-derivative roughness
#'                   penalty. A genuinely smooth reconstruction, outside the
#'                   polygonal family, included as a contrast.
#' @param K dyadic truncation level.
#' @param grid `vertices` evaluates on the 2^K + 1 dyadic points (thesis
#'   -faithful, reproduces the printed tables); `input` evaluates back on the
#'   original sampling grid.
#' @param lambda roughness penalty for `bspline_pen`; NULL selects by GCV over
#'   `RECON_LAMBDA_GCV_GRID`.
#' @return list(X = matrix on the evaluation grid, argvals = that grid,
#'   method, K, grid)
reconstruct <- function(X,
                        method  = RECON_METHOD_DEFAULT,
                        K       = RECON_K_DEFAULT,
                        grid    = RECON_GRID_DEFAULT,
                        lambda  = RECON_LAMBDA_DEFAULT,
                        argvals = NULL) {
  if (!method %in% RECONSTRUCTION_METHODS) {
    stop("unknown reconstruction method: ", method,
         " (expected one of ", paste(RECONSTRUCTION_METHODS, collapse = ", "), ")")
  }
  X <- as.matrix(X)
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = nrow(X))
  out_grid <- recon_grid(method, K, grid, nrow(X))

  if (identical(method, "none")) {
    return(list(X = X, argvals = argvals, method = method, K = NA_integer_,
                grid = "input"))
  }

  Xr <- switch(method,
    fs_interp = {
      # FS truncation == polygonal interpolation at the dyadic vertices.
      V <- apply(X, 2L, function(col) approx(argvals, col, xout = seq(0, 1, by = 2^(-as.integer(K))))$y)
      V <- matrix(V, nrow = 2^as.integer(K) + 1L)
      if (identical(grid, "input")) {
        apply(V, 2L, function(col)
          approx(seq(0, 1, by = 2^(-as.integer(K))), col, xout = argvals)$y)
      } else {
        V
      }
    },
    fs_ls = {
      fd <- Data2fd(argvals, X, basisobj = fs_basis(K))
      eval.fd(out_grid, fd)
    },
    bspline_pen = {
      b <- pen_basis(K)
      lam <- if (is.null(lambda)) .gcv_lambda(argvals, X, b) else lambda
      fdp <- fdPar(b, int2Lfd(2L), lam)
      fd  <- smooth.basis(argvals, X, fdp)$fd
      eval.fd(out_grid, fd)
    })

  Xr <- matrix(as.numeric(Xr), nrow = length(out_grid))
  list(X = Xr, argvals = out_grid, method = method, K = as.integer(K),
       grid = grid)
}

RECON_LAMBDA_GCV_GRID <- 10^seq(-8, 2, by = 1)

#' Pick lambda by pooled GCV over `RECON_LAMBDA_GCV_GRID`.
.gcv_lambda <- function(argvals, X, basis) {
  scores <- vapply(RECON_LAMBDA_GCV_GRID, function(lam) {
    fdp <- fdPar(basis, int2Lfd(2L), lam)
    sum(smooth.basis(argvals, X, fdp)$gcv, na.rm = TRUE)
  }, numeric(1))
  RECON_LAMBDA_GCV_GRID[which.min(scores)]
}

#' Short human-readable label, for CSV columns and plot facets.
recon_label <- function(method, K = NA, grid = NA) {
  if (identical(method, "none")) return("none")
  sprintf("%s_K%d_%s", method, as.integer(K), grid)
}
