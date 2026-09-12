#' m-split aggregation for paper 1's T^(d) scan statistic.
#'
#' Paper 1 (A1) defines the aggregated statistic
#'   T^(d)_{n, rho}(H) = max_{1 <= i <= d} T_{N, rho}(h_i; X_{i, 1}, ..., X_{i, N}),
#' with N = n / d and independent kernels / projections per chunk. For our
#' purposes the "kernel" is a projection that draws a fresh random direction
#' each call (e.g. `project_ito` samples a new Wiener path, `project_c_max_var`
#' samples a new critical point), so each chunk genuinely sees a different
#' feature function.
#'
#' The `m == 1` path short-circuits the split and is equivalent to calling
#' `run(X)` directly.
#'
#' Null-distribution calibration: since the aggregate is max over m
#' (asymptotically iid) single-chunk statistics, use the (1-alpha)^(1/m)
#' upper quantile of the cached null -- `critical_value(null, alpha, m.d = m)`
#' in `shared/weights.R` does this directly.
#'
#' @param X full sample: matrix with n columns (functional) or numeric
#'   length-n vector (scalar).
#' @param m number of splits. Must divide n evenly.
#' @param run closure run(chunk) -> scalar statistic. The closure must
#'   encapsulate projection + statistic, because projection needs to be
#'   re-drawn per chunk (each call generates a fresh random direction where
#'   applicable).
#' @return scalar max statistic over the m chunks.
scan_over_splits <- function(X, m, run) {
  if (m == 1L) return(run(X))
  n <- if (is.matrix(X)) ncol(X) else length(X)
  if (n %% m != 0L) {
    stop(sprintf("scan_over_splits: n = %d not divisible by m = %d", n, m))
  }
  N <- n %/% m
  vals <- numeric(m)
  for (i in seq_len(m)) {
    cols <- ((i - 1L) * N + 1L):(i * N)
    chunk <- if (is.matrix(X)) X[, cols, drop = FALSE] else X[cols]
    vals[i] <- run(chunk)
  }
  max(vals)
}

#' Build a `run(chunk) -> scalar` closure uniform across scalar / matrix
#' projections.
#'
#' For scalar projections (numeric(n) output) the statistic runs on the
#' projected vector directly. For matrix projections (e.g. `project_pc_d`
#' returning matrix(d, n)) the statistic runs on each row and the max is
#' taken -- this is the inner `max_{h in H}` of paper 1's T_{n, rho}(H).
#' The CV must be adjusted separately (see `projection_multiplicity`) to
#' account for the d features in the overall max distribution.
#'
#' `matrix_in = TRUE` is the escape hatch for statistics that consume the
#' full (l + 1) x n functional sample (e.g. Birbilas's E_{n,beta,gamma},
#' which couples the time and t axes). With `matrix_in`, the projection
#' output is forwarded to `stat_fn` wholesale -- no per-row iteration.
make_run <- function(project, stat_fn, gamma, stat_args, matrix_in = FALSE) {
  function(chunk) {
    projected <- project(chunk)
    if (matrix_in) {
      do.call(stat_fn, c(list(projected, gamma), stat_args))
    } else if (is.matrix(projected)) {
      max(vapply(seq_len(nrow(projected)),
                 function(k) do.call(stat_fn,
                                     c(list(projected[k, ], gamma), stat_args)),
                 numeric(1)))
    } else {
      do.call(stat_fn, c(list(projected, gamma), stat_args))
    }
  }
}
