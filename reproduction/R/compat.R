# Compatibility layer between the original pipeline's conventions and the
# fdachange package. The original scripts expect `generate_<SCENARIO>` /
# `cdf_<SCENARIO>` symbols, an `EC_EXACT` global selecting the change-window
# convention, `kernel_fn(name)`, `far1_beta_registry`, and
# `fmt_duration`. This file provides them on top of the package so the
# ported entry points keep their original bodies.
#
# The package generators are bit-identical to the original scenario files
# (proven by the package's golden-value suite, including the EnvStats
# Pareto draw streams), so routing through them preserves every RNG stream.

.ec_exact <- function() isTRUE(get0("EC_EXACT", ifnotfound = FALSE))

generate_S1 <- function(n, a = 0, ec.interval = c(0.4, 0.7)) {
  fdachange::generate_s1(n, a, ec.interval, exact_window = .ec_exact())
}
generate_S2 <- function(n, a = 0, ec.interval = c(0.4, 0.7)) {
  fdachange::generate_s2(n, a, ec.interval, exact_window = .ec_exact())
}
generate_S3 <- function(n, a = 0, ec.interval = c(0.4, 0.7)) {
  fdachange::generate_s3(n, a, ec.interval, exact_window = .ec_exact())
}
generate_S4 <- function(n, a = 0, ec.interval = c(0.4, 0.7)) {
  fdachange::generate_s4(n, a, ec.interval, exact_window = .ec_exact())
}
generate_S5 <- function(n, a = 0, ec.interval = c(0.4, 0.7)) {
  fdachange::generate_s5(n, a, ec.interval, exact_window = .ec_exact())
}
generate_F1 <- function(n, a = 0, ec.interval = c(0.4, 0.7), l = 100) {
  fdachange::generate_f1(n, a, ec.interval, l, exact_window = .ec_exact())
}
generate_F2 <- function(n, a = 0, ec.interval = c(0.4, 0.7), l = 100) {
  fdachange::generate_f2(n, a, ec.interval, l, exact_window = .ec_exact())
}
generate_F3 <- function(n, a = 0, ec.interval = c(0.4, 0.7), l = 100) {
  fdachange::generate_f3(n, a, ec.interval, l, exact_window = .ec_exact())
}
generate_F4 <- function(n, a = 0, ec.interval = c(0.25, 0.45), l = 100,
                        beta_t = fdachange::far1_betas$beta1) {
  fdachange::generate_f4(n, a, ec.interval, l, beta_t)
}
generate_F5 <- function(n, a = 0, ec.interval = c(0.4, 0.7), l = 100,
                        spike_alpha = 0.45) {
  fdachange::generate_f5(n, a, ec.interval, l, spike_alpha,
                         exact_window = .ec_exact())
}

cdf_S1 <- fdachange::scenario_cdf("S1")
cdf_S2 <- fdachange::scenario_cdf("S2")
cdf_S3 <- fdachange::scenario_cdf("S3")
cdf_S4 <- fdachange::scenario_cdf("S4")
cdf_S5 <- fdachange::scenario_cdf("S5")
cdf_F1 <- fdachange::scenario_cdf("F1")
cdf_F2 <- fdachange::scenario_cdf("F2")
cdf_F3 <- fdachange::scenario_cdf("F3")
cdf_F4 <- fdachange::scenario_cdf("F4")
cdf_F5 <- fdachange::scenario_cdf("F5")

kernel_fn <- function(name) {
  switch(name,
    gaussian = fdachange::kernel_gaussian,
    min      = fdachange::kernel_min,
    `2max`   = ,
    `2_max`  = fdachange::kernel_2max,
    stop("unknown kernel: ", name))
}

far1_beta_registry <- fdachange::far1_betas

fmt_duration <- function(s) {
  if (!is.finite(s) || s < 0) return("--")
  if (s < 60)     return(sprintf("%.1fs", s))
  if (s < 3600)   return(sprintf("%dm%02ds", as.integer(s) %/% 60L, as.integer(s) %% 60L))
  if (s < 86400)  return(sprintf("%dh%02dm", as.integer(s) %/% 3600L, (as.integer(s) %% 3600L) %/% 60L))
  sprintf("%dd%02dh", as.integer(s) %/% 86400L, (as.integer(s) %% 86400L) %/% 3600L)
}
