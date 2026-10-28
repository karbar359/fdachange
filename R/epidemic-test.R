#' Epidemic change-point test
#'
#' Runs one of the four scan statistics against its null calibration and
#' returns a standard `htest`-style result with the statistic, a p-value,
#' the level-`alpha` critical value, and the detected change window.
#'
#' The null distribution is resolved automatically from the shipped
#' Monte-Carlo cache for the method ([wilcoxon_null_draws()],
#' [lp_null_draws()], [rkhs_null_draws()], [far1_null_draws()]). Pass
#' your own draws via `null=` to use a scenario-specific or finite-sample
#' calibration instead (required when
#' `kernel` is a custom function, or for RKHS baselines that are not
#' cached).
#'
#' Calibration caveat: the Wilcoxon scan and the ecdf-path Lp scan
#' (`cdf = NULL`) are rank-based, so their nulls apply to any continuous
#' data. The RKHS caches are tied to the H0 distribution named by
#' `baseline`; for data with a different null marginal, simulate the null
#' under your H0 ([simulate_rkhs_null()] with `rdist=`) and pass it via
#' `null=` or [save_null_draws()].
#'
#' With `splits = m > 1` the sample is cut into `m` equal consecutive
#' chunks, the statistic is the maximum over chunks (see
#' [scan_over_splits()]), and the p-value / critical value are adjusted
#' by `F(t)^m`. The detected window is only reported for `splits = 1`.
#'
#' @param x numeric vector for `"wilcoxon"`, `"lp"`, `"rkhs"` (project
#'   functional data first, see [project_l2()]); numeric matrix
#'   `(grid x n)` for `"far1"`.
#' @param method which scan statistic to use.
#' @param gamma Hoelder / window-length weight exponent (scalar).
#' @param beta logarithmic weight exponent (Wilcoxon, RKHS).
#' @param kernel Wilcoxon: `"h_w"` or `"h_0"`; RKHS: kernel name or
#'   function. Default: `"h_w"` for Wilcoxon, `"gaussian"` for RKHS.
#' @param p,cdf,t_grid Lp scan parameters, see [lp_statistic()].
#' @param beta_t Birbilas beta operator, see [far1_statistic()].
#' @param baseline RKHS null baseline tag, see [rkhs_null_draws()].
#' @param alpha significance level for the reported critical value.
#' @param splits number of sample splits (default 1).
#' @param null optional numeric vector of null draws overriding the
#'   automatic calibration.
#' @return an object of class `c("epidemic_htest", "htest")`: a list with
#'   `statistic`, `p.value`, `critical.value`, `alpha`, `interval`
#'   (fractional change window `c(k, m)/n`), `k`, `m`, `parameter`,
#'   `method`, `data.name`, `n.null` (number of null draws backing the
#'   calibration) and `splits`.
#' @references Bartkus (2026); Bartkus and Rackauskas (2025, 2026);
#'   Birbilas (2025).
#' @examples
#' set.seed(1)
#' x <- rnorm(200)
#' x[81:140] <- x[81:140] + 0.8
#' epidemic_test(x, "wilcoxon", gamma = 0)
#' epidemic_test(x, "rkhs", gamma = 0.25, kernel = "gaussian")
#'
#' # functional data through a projection
#' X <- generate_f1(100, a = 3)
#' epidemic_test(project_l2(X), "lp", gamma = 0.25)
#' @export
epidemic_test <- function(x,
                          method = c("wilcoxon", "lp", "rkhs", "far1"),
                          gamma, beta = 0, kernel = NULL, p = 2, cdf = NULL,
                          t_grid = 100, beta_t = far1_betas$beta1,
                          baseline = "", alpha = 0.05, splits = 1L,
                          null = NULL) {
  method <- match.arg(method)
  data.name <- deparse(substitute(x))
  stopifnot(length(gamma) == 1L, splits >= 1L)
  splits <- as.integer(splits)

  scan_fn <- switch(method,
    wilcoxon = function(z) wilcoxon_scan(z, gamma, beta,
                                         kernel = if (is.null(kernel)) "h_w"
                                                  else kernel),
    lp       = function(z) lp_scan(z, gamma, p = p, cdf = cdf,
                                   t_grid = t_grid),
    rkhs     = function(z) rkhs_scan(z, gamma,
                                     kernel = if (is.null(kernel)) "gaussian"
                                              else kernel, beta = beta),
    far1 = function(z) far1_scan(z, gamma, beta_t = beta_t))

  if (splits == 1L) {
    sc <- scan_fn(x)
    stat <- sc$stat; k <- sc$k; m <- sc$m; interval <- sc$interval
  } else {
    stat <- scan_over_splits(x, splits, function(chunk) scan_fn(chunk)$stat)
    k <- NA_integer_; m <- NA_integer_; interval <- c(NA_real_, NA_real_)
  }

  draws <- if (!is.null(null)) {
    stopifnot(is.numeric(null), length(null) > 1)
    null
  } else {
    switch(method,
      wilcoxon = wilcoxon_null_draws(gamma, beta),
      lp       = lp_null_draws(gamma, p),
      rkhs     = {
        if (is.function(kernel)) {
          stop("rkhs with a custom kernel function needs explicit null ",
               "draws: pass null = simulate_rkhs_null(...)", call. = FALSE)
        }
        rkhs_null_draws(gamma, if (is.null(kernel)) "gaussian" else kernel,
                        baseline)
      },
      far1 = far1_null_draws(gamma))
  }
  p.value <- 1 - mean(draws <= stat)^splits
  cv <- critical_value(draws, alpha, splits)
  n.null <- length(draws)

  method_label <- switch(method,
    wilcoxon = "Epidemic change scan test (Wilcoxon-type antisymmetric kernel)",
    lp       = sprintf("Epidemic change scan test (L%s empirical process)", p),
    rkhs     = "Epidemic change scan test (RKHS kernel mean embedding)",
    far1 = "Epidemic change test for FAR(1) (Birbilas E_{n,beta,gamma})")

  parameter <- c(gamma = gamma)
  if (beta != 0) parameter <- c(parameter, beta = beta)
  if (method == "lp") parameter <- c(parameter, p = p)
  if (splits > 1L) parameter <- c(parameter, splits = splits)

  structure(
    list(statistic = c(T = stat), parameter = parameter, p.value = p.value,
         critical.value = cv, alpha = alpha, interval = interval,
         k = k, m = m, splits = splits, n.null = n.null,
         alternative = "an epidemic (changed then restored) segment",
         method = method_label, data.name = data.name),
    class = c("epidemic_htest", "htest"))
}

#' @rdname epidemic_test
#' @param ... further arguments passed to [print()] methods.
#' @export
print.epidemic_htest <- function(x, ...) {
  y <- x
  class(y) <- "htest"
  print(y, ...)
  cat(sprintf("critical value (alpha = %s): %.4f  [%s null draws]\n",
              format(x$alpha), x$critical.value, format(x$n.null)))
  if (!is.na(x$k)) {
    cat(sprintf("detected window: observations %d..%d (fractions %.3f..%.3f)\n",
                x$k + 1L, x$m, x$interval[1], x$interval[2]))
  }
  invisible(x)
}
