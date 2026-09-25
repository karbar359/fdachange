# The package must reproduce the audited original implementation to
# machine precision: statistics, scan windows, critical values from every
# shipped cache cell, scenario RNG streams, and projections.

GAMMAS <- c(0, 0.25, 0.45)
NS <- c(40L, 60L, 120L)
SEED <- 7L

test_that("statistics match the golden capture", {
  for (n in NS) {
    set.seed(SEED + n); x <- rnorm(n)
    for (g in GAMMAS) {
      for (kern in c("h_w", "h_0")) {
        expect_golden(sprintf("wilcoxon|n=%d|g=%s|k=%s", n, g, kern),
                      wilcoxon_statistic(x, g, kernel = kern))
      }
      for (p in c(2L, 4L)) {
        expect_golden(sprintf("lp|n=%d|g=%s|p=%d|cdf=pnorm", n, g, p),
                      lp_statistic(x, g, p = p, cdf = pnorm))
        expect_golden(sprintf("lp|n=%d|g=%s|p=%d|cdf=ecdf", n, g, p),
                      lp_statistic(x, g, p = p))
      }
      for (kern in c("gaussian", "min", "2max")) {
        expect_golden(sprintf("rkhs|n=%d|g=%s|k=%s", n, g, kern),
                      rkhs_statistic(x, g, kernel = kern))
      }
    }
  }
})

test_that("Birbilas statistic matches the golden capture", {
  for (n in c(40L, 60L)) {
    set.seed(SEED + n); Xf <- generate_f4(n, 0, c(0.4, 0.6))
    for (g in c(0.25, 0.4)) {
      for (b in c("beta1", "beta3")) {
        expect_golden(sprintf("birbilas|n=%d|g=%s|b=%s|a=0", n, g, b),
                      birbilas_statistic(Xf, g, beta_t = birbilas_betas[[b]]))
      }
    }
    set.seed(SEED + n); Xf2 <- generate_f4(n, 2, c(0.25, 0.45))
    expect_golden(sprintf("birbilas|n=%d|g=0.4|b=beta1|a=2", n),
                  birbilas_statistic(Xf2, 0.4))
  }
})

test_that("scan windows match the golden capture", {
  set.seed(SEED + 60L); x60 <- rnorm(60)
  sc <- wilcoxon_scan(x60, 0.25)
  expect_golden("scan|wilcoxon|stat", sc$stat)
  expect_golden("scan|wilcoxon|k", sc$k)
  expect_golden("scan|wilcoxon|m", sc$m)
  sc <- lp_scan(x60, 0.25, cdf = pnorm)
  expect_golden("scan|lp|stat", sc$stat)
  expect_golden("scan|lp|k", sc$k)
  expect_golden("scan|lp|m", sc$m)
  sc <- rkhs_scan(x60, 0.25, kernel = "gaussian")
  expect_golden("scan|rkhs|stat", sc$stat)
  expect_golden("scan|rkhs|k", sc$k)
  expect_golden("scan|rkhs|m", sc$m)
})

test_that("critical values from every shipped cache cell match", {
  w_cells <- rbind(
    expand.grid(g = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                      0.495), b = 0),
    data.frame(g = 0.5, b = 1))
  for (i in seq_len(nrow(w_cells))) {
    g <- w_cells$g[i]; b <- w_cells$b[i]
    for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
      expect_golden(sprintf("cv|wilcoxon|g=%s|b=%s|a=%s|s=%d", g, b, al, m),
                    wilcoxon_critical_value(g, b, alpha = al, splits = m))
    }
  }
  for (g in GAMMAS) for (p in c(2L, 4L, 10L)) {
    for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
      expect_golden(sprintf("cv|lp|g=%s|p=%d|a=%s|s=%d", g, p, al, m),
                    lp_critical_value(g, p, alpha = al, splits = m))
    }
  }
  rk_cells <- list(c("gaussian", ""), c("min", "rnorm.sq"), c("min", "runif"),
                   c("2max", "runif"), c("2max", "rbeta.1.9"))
  for (g in GAMMAS) for (cell in rk_cells) {
    for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
      expect_golden(sprintf("cv|rkhs|g=%s|k=%s|b=%s|a=%s|s=%d",
                            g, cell[1], cell[2], al, m),
                    rkhs_critical_value(g, cell[1], cell[2],
                                        alpha = al, splits = m))
    }
  }
  for (g in c(0, 0.25, 0.4, 0.45)) {
    for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
      expect_golden(sprintf("cv|birbilas|g=%s|a=%s|s=%d", g, al, m),
                    birbilas_critical_value(g, alpha = al, splits = m))
    }
  }
})

test_that("scenario RNG streams match the original (EnvStats parity)", {
  gens <- list(s1 = generate_s1, s2 = generate_s2, s3 = generate_s3,
               s4 = generate_s4, s5 = generate_s5)
  for (s in names(gens)) {
    set.seed(11L)
    expect_golden(sprintf("stream|%s|a=1|exact", s),
                  sum(gens[[s]](50, a = 1, exact_window = TRUE)))
  }
  set.seed(11L)
  expect_golden("stream|f1|a=3|exact", sum(generate_f1(20, a = 3, l = 20)))
  set.seed(11L)
  expect_golden("stream|f3|a=3|exact", sum(generate_f3(20, a = 3, l = 20)))
  set.seed(11L)
  expect_golden("stream|f2|a=0.5|exact", sum(generate_f2(20, a = 0.5, l = 20)))
  set.seed(11L)
  expect_golden("stream|f4|a=2", sum(generate_f4(20, a = 2, l = 20)))
  set.seed(11L)
  expect_golden("stream|s1|a=1|legacy",
                sum(generate_s1(50, a = 1, exact_window = FALSE)))
  set.seed(11L)
  expect_golden("stream|s3|a=1|legacy",
                sum(generate_s3(50, a = 1, exact_window = FALSE)))
})

test_that("projections match the golden capture", {
  set.seed(11L); Xp <- generate_f1(20, a = 0, l = 20)
  expect_golden("proj|l2|sum", sum(project_l2(Xp)))
  set.seed(13L)
  expect_golden("proj|ito|sum", sum(project_ito(Xp)))
  skip_if_not_installed("fda")
  skip_if_not_installed("fda.usc")
  set.seed(13L)
  expect_golden("proj|c_max_var|sum", sum(project_c_max_var(Xp)),
                tolerance = 1e-8)
  set.seed(13L)
  expect_golden("proj|wiener|sum", sum(project_wiener(Xp)), tolerance = 1e-8)
  set.seed(13L)
  # fda.usc's fdata2pc warns about a missing %dopar% backend; harmless here
  expect_golden("proj|pc1|sumabs", sum(abs(suppressWarnings(project_pc1(Xp)))),
                tolerance = 1e-8)
})
