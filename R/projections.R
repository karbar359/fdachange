# Functional-to-scalar projections. project_l2() and project_ito() are
# base R; the remaining projections expand the curves in a B-spline basis
# via the fda / fda.usc packages (Suggests) and error informatively when
# those are not installed.

.check_fda <- function(need_usc = TRUE) {
  ok <- requireNamespace("fda", quietly = TRUE) &&
    (!need_usc || requireNamespace("fda.usc", quietly = TRUE))
  if (!ok) {
    stop("this projection needs the ",
         if (need_usc) "'fda' and 'fda.usc' packages" else "'fda' package",
         "; install with install.packages(",
         if (need_usc) "c(\"fda\", \"fda.usc\")" else "\"fda\"", ")",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Order-2 B-spline basis with dyadic knots at J = 4 (17 knots on [0, 1]):
# the Schauder basis of Bartkus and Rackauskas (2025). Cached on first use.
.basis_env <- new.env(parent = emptyenv())
.default_basis <- function() {
  if (is.null(.basis_env$b)) {
    knots <- seq(0, 1, by = 2^(-4L))
    .basis_env$b <- fda::create.bspline.basis(norder = 2L, breaks = knots)
  }
  .basis_env$b
}

# Wrap a sample matrix X (grid rows, one curve per column) as an fda.usc
# fdata object, expanding in a B-spline basis over [0, 1].
.matrix_to_fdata <- function(X, argvals = NULL, basis = NULL) {
  .check_fda()
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = nrow(X))
  if (is.null(basis))   basis   <- .default_basis()
  X.fd <- fda::Data2fd(argvals, X, basisobj = basis)
  fda.usc::fdata(X.fd)
}

#' Functional-to-scalar projections
#'
#' Reduce a matrix of sampled curves (rows = evaluation grid on `[0, 1]`,
#' one observation per **column**) to a scalar series of length `n` that
#' the scan statistics consume:
#'
#' * `project_l2()`: the L2 norm \eqn{(\int_0^1 X_j(t)^2 dt)^{1/2}} by the
#'   trapezoidal rule. Deterministic per curve; base R.
#' * `project_ito()`: midpoint (Stratonovich-type) integral of every curve
#'   against one shared random Wiener path per call; base R.
#' * `project_c_max_var()`: evaluates every curve at a randomly chosen
#'   critical point of the functional variance (the projection of Bartkus
#'   and Rackauskas, 2025). Needs fda + fda.usc.
#' * `project_wiener()`: L2 inner product with one random Brownian path.
#'   Needs fda + fda.usc.
#' * `project_pc1()`: inner product with the leading functional principal
#'   component. Needs fda + fda.usc.
#' * `project_pc_d()`: inner products with the leading `d` components;
#'   returns a `d x n` matrix (one scalar series per row) for the
#'   max-aggregation of [make_run()], with critical values adjusted by
#'   [projection_multiplicity()]. Needs fda + fda.usc.
#'
#' The stochastic projections (`ito`, `c_max_var`, `wiener`) draw a fresh
#' random direction on every call; that is what [scan_over_splits()]
#' relies on. Projections consume whatever grid they are handed and do
#' not smooth: reconstruct first with [reconstruct()] when the analysis
#' calls for it.
#'
#' @param X numeric matrix of curves, `(grid points) x n`.
#' @param argvals evaluation grid (default: equispaced on `[0, 1]`).
#' @param basis fda basis object for the fda-backed projections (default:
#'   order-2 B-splines on dyadic knots, `J = 4`).
#' @param d number of principal components for `project_pc_d()`.
#' @param ... unused; accepted so all projections share one signature.
#' @return numeric vector of length `n` (`project_pc_d()`: a `d x n`
#'   matrix).
#' @examples
#' set.seed(1)
#' X <- generate_f1(50, a = 0)
#' project_l2(X)[1:5]
#' project_ito(X)[1:5]
#' @export
project_l2 <- function(X, ...) {
  L <- nrow(X)
  sq <- X^2
  w <- rep(1, L); w[1L] <- w[L] <- 0.5
  integ <- (w %*% sq) / (L - 1L)
  sqrt(as.numeric(integ))
}

#' @rdname project_l2
#' @export
project_ito <- function(X, ...) {
  L <- nrow(X)
  l <- L - 1L
  X_shift <- rbind(0, X[-L, , drop = FALSE])
  B_diffs <- rnorm(L, 0, 1 / l)
  as.numeric(B_diffs %*% ((X + X_shift) / 2))
}

#' @rdname project_l2
#' @export
project_c_max_var <- function(X, argvals = NULL, basis = NULL) {
  X.fdata <- .matrix_to_fdata(X, argvals, basis)
  projection_point <- runif(1, 0, 1)
  d <- fda.usc::fdata.deriv(fda.usc::func.var(X.fdata), method = "diff")$data
  t.grid <- X.fdata$argvals
  idx <- which(d[-length(d)] * d[-1] < 0)
  if (length(idx) > 0L) {
    projection_point <- t.grid[idx[ceiling(projection_point * length(idx))]]
  }
  as.numeric(fda::eval.fd(fdobj = fda.usc::fdata2fd(X.fdata),
                          evalarg = projection_point))
}

#' @rdname project_l2
#' @export
project_wiener <- function(X, argvals = NULL, basis = NULL) {
  X.fdata <- .matrix_to_fdata(X, argvals, basis)
  t.grid  <- X.fdata$argvals
  W.fdata <- fda.usc::rproc2fdata(1L, t = t.grid, sigma = "wiener")
  as.numeric(fda.usc::inprod.fdata(X.fdata, W.fdata))
}

#' @rdname project_l2
#' @export
project_pc1 <- function(X, argvals = NULL, basis = NULL) {
  X.fdata <- .matrix_to_fdata(X, argvals, basis)
  pca <- fda.usc::fdata2pc(X.fdata, ncomp = 1L)
  as.numeric(fda.usc::inprod.fdata(X.fdata, pca$rotation[1]))
}

#' @rdname project_l2
#' @export
project_pc_d <- function(X, d = 3L, argvals = NULL, basis = NULL) {
  X.fdata <- .matrix_to_fdata(X, argvals, basis)
  pca <- fda.usc::fdata2pc(X.fdata, ncomp = d)
  res <- matrix(NA_real_, nrow = d, ncol = ncol(X))
  for (k in seq_len(d)) {
    res[k, ] <- as.numeric(fda.usc::inprod.fdata(X.fdata, pca$rotation[k]))
  }
  res
}

#' Resolve a projection by name
#'
#' Maps a projection id to its function: `"none"` (identity, for scalar
#' data), `"l2"`, `"ito"`, `"c_max_var"`, `"wiener"`, `"pc1"`, `"pc_d"`.
#'
#' @param name projection id.
#' @param pc_d_n number of components for the `"pc_d"` projection, used
#'   by `projection_multiplicity()`.
#' @return `projection_fn()`: the projection function.
#' @examples
#' proj <- projection_fn("l2")
#' set.seed(1)
#' proj(generate_f1(30, a = 0))[1:3]
#' @export
projection_fn <- function(name) {
  switch(name,
    none      = function(X, ...) X,
    c_max_var = project_c_max_var,
    ito       = project_ito,
    l2        = project_l2,
    wiener    = project_wiener,
    pc1       = project_pc1,
    pc_d      = project_pc_d,
    stop("unknown projection: ", name))
}

#' @rdname projection_fn
#' @return `projection_multiplicity()`: the number of scalar series the
#'   projection emits (adjust critical values by
#'   `critical_value(draws, alpha, splits * multiplicity)`).
#' @export
projection_multiplicity <- function(name, pc_d_n = 3L) {
  switch(name,
    pc_d = as.integer(pc_d_n),
    1L)
}
