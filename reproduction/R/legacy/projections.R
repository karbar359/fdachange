#' Functional projections.
#'
#' Reduce a matrix(l+1, n) of sampled curves to a numeric(n) scalar series
#' that the existing scalar scan statistics (RKHS / Lp / Wilcoxon) can
#' consume unchanged. Scalar scenarios use `none` (identity); functional
#' scenarios pass their `matrix` output through one of the named
#' projections below.
#'
#' Paper 1 (A1) uses `c_max_var`. Other projections (Ito, Wiener, PC1/PC2,
#' PC_d, L2, Normally_weighted_mean) are dormant in the legacy code and
#' can be added one-by-one as needed.
#'
#' SMOOTHING IS THE CALLER'S JOB. These functions consume whatever grid they
#' are handed; they do not reconstruct. That matters because they are not
#' uniform about it internally: `project_c_max_var`, `project_wiener`,
#' `project_pc1` and `project_pc_d` wrap their input in an fd object via
#' `matrix_to_fdata` (and so implicitly least-squares it onto the K = 4
#' polygonal basis), while `project_l2` and `project_ito` work directly on the
#' sampled grid. The dissertation's application (printed Tables 7.1 and 7.2)
#' was computed on RECONSTRUCTED curves for all three reported projections, so
#' passing a raw matrix here does not reproduce it: at gamma = 0.45 the raw
#' L2 path gives T = 2.17 against a printed 3.44.
#'
#' Reconstruct first with `shared/reconstruction.R::reconstruct()`, then
#' project. See `applications/litgrid/reproduce_tables.R` for the exact
#' combination that reproduces the printed tables, and
#' `smoothing_sensitivity.R` for how much the choice matters.

suppressPackageStartupMessages({
  library(fda)
  library(fda.usc)
})

# Paper 1's Schauder B-spline basis: order 2, dyadic knots at J = 4 (17 knots
# on [0, 1]). Cached on first use because building the basis is cheap but
# repetitive across replicates.
.default_basis <- local({
  b <- NULL
  function() {
    if (is.null(b)) {
      knots <- seq(0, 1, by = 2^(-4L))
      b <<- create.bspline.basis(norder = 2L, breaks = knots)
    }
    b
  }
})

#' Wrap a sample matrix X (l+1 rows, one curve per column) as an `fdata`
#' object, expanding in a B-spline basis over [0, 1].
matrix_to_fdata <- function(X, argvals = NULL, basis = NULL) {
  if (is.null(argvals)) argvals <- seq(0, 1, length.out = nrow(X))
  if (is.null(basis))   basis   <- .default_basis()
  X.fd <- Data2fd(argvals, X, basisobj = basis)
  fda.usc::fdata(X.fd)
}

#' C_max_var projection (paper 1 A1).
#'
#' Computes the functional variance, differentiates it, finds sign changes
#' of the derivative (critical points of the variance), and picks one
#' uniformly at random. Every curve is then evaluated at that point,
#' yielding a scalar series of length n.
#'
#' If no critical point exists (rare -- constant derivative sign), falls
#' back to a uniform random point in (0, 1).
project_c_max_var <- function(X, argvals = NULL, basis = NULL) {
  X.fdata <- matrix_to_fdata(X, argvals, basis)
  projection_point <- runif(1, 0, 1)
  d <- fdata.deriv(func.var(X.fdata), method = "diff")$data
  t.grid <- X.fdata$argvals
  idx <- which(d[-length(d)] * d[-1] < 0)
  if (length(idx) > 0L) {
    projection_point <- t.grid[idx[ceiling(projection_point * length(idx))]]
  }
  as.numeric(eval.fd(fdobj = fdata2fd(X.fdata), evalarg = projection_point))
}

#' Ito (midpoint Stratonovich) projection.
#'
#' Integrates every curve against a single random Wiener path:
#'   X_projected_j = sum_i B_i * (X[i, j] + X[i-1, j]) / 2
#' with B_i ~ N(0, 1/l^2) -- matches the legacy R implementation in
#' `helper_functions.R:368` (`B.diffs = rnorm(l+1, 0, 1/l)`). One path is
#' drawn per call, shared across all n curves in the sample. Used in
#' paper 1 tables 1-3 ("Ito integral" rows).
project_ito <- function(X, ...) {
  L <- nrow(X)
  l <- L - 1L
  X_shift <- rbind(0, X[-L, , drop = FALSE])
  B_diffs <- rnorm(L, 0, 1 / l)
  as.numeric(B_diffs %*% ((X + X_shift) / 2))
}

#' L^2 norm projection.
#'
#' X_projected_j = (int_0^1 X(t, j)^2 dt)^{1/2}, computed by the trapezoidal
#' rule on the grid. Used in paper 1 table 1 ("||.||_L2" rows). This is a
#' feature, not a function of a random direction -- deterministic per curve.
project_l2 <- function(X, ...) {
  L <- nrow(X)
  sq <- X^2
  w <- rep(1, L); w[1L] <- w[L] <- 0.5
  integ <- (w %*% sq) / (L - 1L)
  sqrt(as.numeric(integ))
}

#' Wiener (Brownian) projection.
#'
#' Draws one realisation W of standard Brownian motion on the argvals grid
#' (via fda.usc::rproc2fdata with sigma = "wiener"), takes the L2 inner
#' product <X_j, W> for every curve in the sample, yields numeric(n). One
#' path per call -- same stochastic-projection discipline as `project_ito`
#' and `project_c_max_var`. Matches `helper_functions.R:23` trajectory
#' generator + line 404 ("Wiener" branch in the PROJECTION.TYPE.mapings).
project_wiener <- function(X, argvals = NULL, basis = NULL) {
  X.fdata <- matrix_to_fdata(X, argvals, basis)
  t.grid  <- X.fdata$argvals
  W.fdata <- rproc2fdata(1L, t = t.grid, sigma = "wiener")
  as.numeric(inprod.fdata(X.fdata, W.fdata))
}

#' PC1 projection: functional-PCA first component.
#'
#' Computes functional PCA on the sample (fda.usc::fdata2pc), takes the
#' leading eigen-function PC_1, and returns the inner product <X_j, PC_1>
#' per curve. Scalar output: numeric(n). Matches paper 1's dormant "PC1"
#' branch (`helper_functions.R:416`).
project_pc1 <- function(X, argvals = NULL, basis = NULL) {
  X.fdata <- matrix_to_fdata(X, argvals, basis)
  pca <- fdata2pc(X.fdata, ncomp = 1L)
  as.numeric(inprod.fdata(X.fdata, pca$rotation[1]))
}

#' PC_d projection: inner product with the leading `d` PCs.
#'
#' Returns a `matrix(d, n)` whose row k holds the scalar series
#' <X_j, PC_k>. Feeds the outer-max aggregation in paper 1's
#' T_{n, rho}(H) = max_{h in H} T(h) with H = {kernel composed with PC_k
#' feature, k = 1..d}. Matches `helper_functions.R:460-470` ("PC_d"
#' branch). Default d = 3 is what paper 1's commented block uses.
#'
#' The caller is responsible for taking the max over rows AND adjusting
#' the critical value by the multiplicity d (see driver / `make_run`).
project_pc_d <- function(X, d = 3L, argvals = NULL, basis = NULL) {
  X.fdata <- matrix_to_fdata(X, argvals, basis)
  pca <- fdata2pc(X.fdata, ncomp = d)
  res <- matrix(NA_real_, nrow = d, ncol = ncol(X))
  for (k in seq_len(d)) {
    res[k, ] <- as.numeric(inprod.fdata(X.fdata, pca$rotation[k]))
  }
  res
}

#' Resolve a projection name to its function.
#' `none` is an identity -- used by scalar scenarios so the driver's
#' project(X) call is a no-op. `pc_d` has a `d` argument; the driver
#' wraps it with PC_D_N.
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

#' Number of underlying scalar series a projection emits.
#' Used to adjust the critical value by (1-alpha)^(1/(M * mult)).
projection_multiplicity <- function(name, pc_d_n = 3L) {
  switch(name,
    pc_d = as.integer(pc_d_n),
    1L)
}
