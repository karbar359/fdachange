# Condensed method x scenario smoke: empirical size in a wide Monte-Carlo
# band at a = 0 and near-full power at a strong alternative. Slow;
# run with FDACHANGE_RUN_SLOW=1 (CI runs it on a schedule, not per push).

skip_if(!nzchar(Sys.getenv("FDACHANGE_RUN_SLOW")),
        "slow: set FDACHANGE_RUN_SLOW=1")

NREP <- 40L
N <- 200L

reject_rate <- function(gen, method, gamma, a, ...) {
  mean(vapply(seq_len(NREP), function(i) {
    x <- gen(N, a = a)
    ht <- epidemic_test(x, method, gamma = gamma, ...)
    ht$p.value < 0.05
  }, logical(1)))
}

test_that("wilcoxon on S1: level holds, strong shift detected", {
  set.seed(101)
  expect_lte(reject_rate(generate_s1, "wilcoxon", 0.25, a = 0), 0.2)
  set.seed(102)
  expect_gte(reject_rate(generate_s1, "wilcoxon", 0.25, a = 2), 0.9)
})

test_that("lp on S1 (parametric CDF): level holds, strong shift detected", {
  set.seed(103)
  expect_lte(reject_rate(generate_s1, "lp", 0.25, a = 0,
                         cdf = scenario_cdf("S1")), 0.2)
  set.seed(104)
  expect_gte(reject_rate(generate_s1, "lp", 0.25, a = 2,
                         cdf = scenario_cdf("S1")), 0.9)
})

test_that("birbilas on F4: level holds, strong jump detected", {
  set.seed(105)
  gen <- function(n, a) generate_f4(n, a = a, l = 50)
  expect_lte(reject_rate(gen, "birbilas", 0.4, a = 0), 0.2)
  set.seed(106)
  expect_gte(reject_rate(gen, "birbilas", 0.4, a = 2.5), 0.9)
})
