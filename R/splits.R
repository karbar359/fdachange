#' Split-sample aggregation of a scan statistic
#'
#' The aggregated statistic of Bartkus and Rackauskas (2025):
#' \deqn{T^{(m)} = \max_{1 \le i \le m} T(h_i; X_{(i-1)N+1}, \dots, X_{iN}),
#'   \qquad N = n / m,}
#' the maximum of the statistic over `m` equal consecutive chunks, with an
#' independently re-drawn random feature (projection) per chunk where the
#' projection is stochastic. Since the aggregate is a maximum over `m`
#' asymptotically independent copies, calibrate with the
#' `(1 - alpha)^(1/m)` upper quantile of the single-chunk null:
#' `critical_value(draws, alpha, splits = m)`.
#'
#' @param X full sample: numeric vector, or matrix with one observation
#'   per column (functional).
#' @param m number of splits; must divide the number of observations.
#' @param run closure `run(chunk) -> scalar statistic`. Encapsulate
#'   projection + statistic in it, so a stochastic projection is re-drawn
#'   per chunk.
#' @return scalar: the maximum statistic over the `m` chunks.
#' @seealso [make_run()], [critical_value()], [epidemic_test()]
#' @examples
#' set.seed(1)
#' X <- generate_f1(120, a = 3)
#' run <- make_run(project_ito, wilcoxon_statistic, gamma = 0.25)
#' scan_over_splits(X, 2, run)
#' @export
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

#' Build a projection + statistic closure for split scans
#'
#' Returns a `run(chunk) -> scalar` closure uniform across scalar and
#' matrix-valued projections, for use with [scan_over_splits()].
#'
#' For scalar projections (numeric output) the statistic runs on the
#' projected vector. For matrix projections (e.g. [project_pc_d()],
#' returning one scalar series per row) the statistic runs on each row and
#' the maximum is taken; adjust the critical value by the multiplicity
#' (see [projection_multiplicity()]). `matrix_in = TRUE` forwards the
#' projection output to `stat_fn` wholesale, for statistics that consume
#' the full functional sample (e.g. [birbilas_statistic()]).
#'
#' @param project projection function, e.g. [project_ito()]; use
#'   `identity` for scalar data.
#' @param stat_fn statistic called as `stat_fn(projected, gamma, ...)`.
#' @param gamma weight exponent forwarded to `stat_fn`.
#' @param stat_args list of further arguments to `stat_fn`.
#' @param matrix_in forward the projection output as one matrix argument.
#' @return function `run(chunk) -> scalar`.
#' @examples
#' run <- make_run(project_l2, rkhs_statistic, gamma = 0.25,
#'                 stat_args = list(kernel = "gaussian"))
#' set.seed(1)
#' run(generate_f1(60, a = 0))
#' @export
make_run <- function(project, stat_fn, gamma, stat_args = list(),
                     matrix_in = FALSE) {
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
