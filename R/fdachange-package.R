#' @keywords internal
#' @details
#' fdachange implements weighted scan tests for the epidemic
#' (changed-segment) alternative: a sample \eqn{X_1, \dots, X_n} is
#' i.i.d. except possibly on an unknown window \eqn{I^* = \{k^*+1, \dots,
#' k^*+l^*\}} where the distribution changes and then changes back. Four
#' test families are provided, each with its own null calibration:
#'
#' * [wilcoxon_statistic()] / [wilcoxon_scan()]: antisymmetric-kernel
#'   (Wilcoxon-type) scan with symmetric Hoelder weighting.
#' * [lp_statistic()] / [lp_scan()]: Lp-norm scan of the weighted
#'   empirical process.
#' * [rkhs_statistic()] / [rkhs_scan()]: kernel mean-embedding (RKHS)
#'   scan.
#' * [birbilas_statistic()] / [birbilas_scan()]: the CUSUM-type
#'   functional statistic of Birbilas (2025).
#'
#' [epidemic_test()] wraps any of the four into a standard `htest`-style
#' object with a p-value, a critical value, and the detected window.
#' Critical values come from shipped Monte-Carlo caches
#' ([wilcoxon_null_draws()] and friends), from the closed-form Kuiper law
#' where it applies ([kuiper_cdf()]), or from user-run simulators
#' ([simulate_wilcoxon_null()] and friends).
#'
#' Functional (curve-per-observation) data enter either through the
#' Birbilas statistic directly, or through a functional-to-scalar
#' projection ([project_l2()], [project_ito()], [project_c_max_var()],
#' ...) optionally preceded by curve reconstruction ([reconstruct()]).
#' The benchmark data-generating scenarios of the dissertation are
#' available as [generate_s1()] ... [generate_f4()].
#'
#' @references
#' Bartkus, K. (2026). *Analysis of Epidemic Changes in Functional Data*.
#' PhD dissertation, Vilnius University.
#'
#' Bartkus, K. and Rackauskas, A. (2025). Lithuanian Mathematical
#' Journal.
#'
#' Bartkus, K. and Rackauskas, A. (2026). Nonlinear Analysis: Modelling
#' and Control.
#'
#' Birbilas, A. (2025). *Functional Data Analysis in Telecommunications*.
#' PhD dissertation, Vilnius University. \doi{10.15388/vu.thesis.841}
#'
#' @importFrom stats approx ecdf pbeta pchisq pnorm quantile rbeta rnorm
#'   runif uniroot
#' @importFrom utils flush.console
"_PACKAGE"
