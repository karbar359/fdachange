#' Simulated FAR(1) telco usage curves
#'
#' Daily voice-usage curves for 11 simulated customer clusters, 50 months
#' each, following a functional autoregressive FAR(1) process. Clusters
#' `cl2`-`cl4` contain a designed epidemic change; `cl5`-`cl11` have no
#' change but trends, level shifts or outliers, for checking how the tests
#' behave when \eqn{H_0} is violated in other ways.
#'
#' Each cluster is simulated as
#' \deqn{Y_j(t) = f_j \big(\mu(t) + Z_j(t)\big) + c_j\,\bar\mu\,
#'   \frac{\sqrt{t}}{\overline{\sqrt{t}}}, \qquad
#'   Z_j(t) = \beta(t) Z_{j-1}(t) + \sqrt{1 - \beta(t)^2}\,
#'   \big(a_j + e_j(t)\big),}
#' clipped at 0, with \eqn{\beta(t) = 0.7\sqrt{t}} (`far1_betas$beta1`),
#' \eqn{a_j \sim N(0, \sigma_a^2)} a month-level and
#' \eqn{e_j(t) \sim N(0, \sigma_e^2)} a day-level innovation, and 100
#' burn-in months. The daily profile \eqn{\mu} and \eqn{\sigma_a,
#' \sigma_e} were fitted on cluster-level averages of a telco
#' bundle-consumption sample; the \eqn{\sqrt{1 - \beta^2}} factor keeps
#' each day's variance at the fitted level, so the FAR(1) process only adds
#' correlation between months (lag 1: 0 on day 1 up to 0.7 on day 28).
#' \eqn{f_j} is a level multiplier and \eqn{c_j} the size of a
#' \eqn{\sqrt{t}}-shaped jump (as in [generate_f4()]), both stored as
#' attributes. The build script is `data-raw/telco_far1.R` in the source
#' repository.
#'
#' | Cluster | Design |
#' |---|---|
#' | `cl1` | no change |
#' | `cl2` | \eqn{\sqrt{t}} jump, +25% of the level (month average), months 37-46 |
#' | `cl3` | \eqn{\sqrt{t}} jump, -25%, months 37-49; outlier months 37, 45, 49 at 1.5 x; 3 outlier days at 3 x |
#' | `cl4` | \eqn{\sqrt{t}} jump from month 37, deepening 2.5% per month |
#' | `cl5` | no jump; rising trend, level 1.00 to 1.30 |
#' | `cl6` | no jump; falling trend, 1.00 to 0.70 |
#' | `cl7` | no jump; higher level, 1.50 |
#' | `cl8` | no jump; lower level, 0.60 |
#' | `cl9` | no jump; 64 outlier days (5%, 2-4 x) and outlier months 12, 16, 30, 32 at 1.5 x |
#' | `cl10` | no jump; high and falling, 1.50 to 1.05 |
#' | `cl11` | no jump; low and rising, 0.60 to 0.84 |
#'
#' `cl1`-`cl4` each use their own fitted model; `cl5`-`cl11` use `cl2`'s.
#' `cl1` was fitted on only four distinct curves and varies much less than
#' the others.
#'
#' Calibration caveat: months are serially dependent, so the shipped null
#' caches (which assume independent observations) and permutation nulls are
#' not exact for this data. Simulate the null from the same FAR(1) model
#' and pass it via `null=` in [epidemic_test()].
#'
#' @format A named list of 11 numeric matrices `cl1`, ..., `cl11`, each
#'   50 x 28: rows are consecutive months, columns are days 1-28 of the
#'   monthly bundle period (the shortest period length), values are average
#'   daily voice minutes. Transpose (`t()`) to get the package's
#'   `(grid x n)` layout. Each matrix carries the attributes
#'   \describe{
#'     \item{`change_months`}{integer months of the designed change window
#'       (empty for no change).}
#'     \item{`multiplier`}{numeric length 50, the level multiplier
#'       \eqn{f_j}.}
#'     \item{`jump`}{numeric length 50, the \eqn{\sqrt{t}} jump size
#'       \eqn{c_j} as a fraction of the level (month average).}
#'     \item{`outliers`}{integer matrix with columns `month`, `day` of the
#'       outlier days.}
#'     \item{`outlier_months`}{integer months with a raised whole-month
#'       level.}
#'   }
#' @source Simulated; model parameters fitted on an anonymised telco
#'   bundle-consumption sample (cluster averages only).
#' @examples
#' names(telco_far1)
#' dim(telco_far1$cl2)
#' attr(telco_far1$cl2, "change_months")
#'
#' # monthly mean usage: the jump in months 37-46
#' plot(rowMeans(telco_far1$cl2), type = "l", xlab = "month",
#'      ylab = "mean daily minutes")
#'
#' # curves in columns, one scalar per month, then scan
#' X <- t(telco_far1$cl2)
#' epidemic_test(project_l2(X), "wilcoxon", gamma = 0.25)
"telco_far1"
