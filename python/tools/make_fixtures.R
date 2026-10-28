#!/usr/bin/env Rscript
# Generate the cross-language parity fixtures from the INSTALLED fdachange
# R package (the reference implementation): fixed input data plus expected
# statistics, scan windows, critical values, and reconstruction outputs,
# written at full double precision (%.17g).
#
# Run from anywhere:  Rscript python/tools/make_fixtures.R

suppressPackageStartupMessages(library(fdachange))

out_dir <- file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "tests", "parity", "fixtures")
out_dir <- normalizePath(out_dir, mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

g17 <- function(x) formatC(x, format = "g", digits = 17)
write_vec <- function(x, name) {
  writeLines(g17(x), file.path(out_dir, name))
}
write_mat <- function(m, name) {
  write.table(matrix(g17(m), nrow = nrow(m)), file.path(out_dir, name),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

rows <- list()
put <- function(id, value) {
  rows[[length(rows) + 1L]] <<- data.frame(id = id, value = g17(value))
}

## ---- fixed inputs -----------------------------------------------------------
set.seed(7);  x60 <- rnorm(60)
set.seed(7);  F4  <- generate_f4(40, 0, c(0.4, 0.6))
set.seed(9);  F5  <- generate_f5(50, a = 3)
write_vec(x60, "input_x60.csv")
write_mat(F4,  "input_F4.csv")
write_mat(F5,  "input_F5.csv")

## ---- statistics and scans ---------------------------------------------------
for (g in c(0, 0.25, 0.45)) {
  put(sprintf("wilcoxon|g=%s", g), wilcoxon_statistic(x60, g))
  put(sprintf("lp_pnorm|g=%s", g), lp_statistic(x60, g, p = 2, cdf = pnorm))
  put(sprintf("lp_ecdf|g=%s", g),  lp_statistic(x60, g, p = 2))
  put(sprintf("rkhs_gaussian|g=%s", g),
      rkhs_statistic(x60, g, kernel = "gaussian"))
}
put("wilcoxon_h0|g=0.25", wilcoxon_statistic(x60, 0.25, kernel = "h_0"))
put("lp_p4|g=0.25", lp_statistic(x60, 0.25, p = 4, cdf = pnorm))
put("rkhs_min|g=0.25",  rkhs_statistic(x60, 0.25, kernel = "min"))
put("rkhs_2max|g=0.25", rkhs_statistic(x60, 0.25, kernel = "2max"))

sc <- wilcoxon_scan(x60, 0.25)
put("scan_wilcoxon|stat", sc$stat); put("scan_wilcoxon|k", sc$k)
put("scan_wilcoxon|m", sc$m)
sc <- lp_scan(x60, 0.25, cdf = pnorm)
put("scan_lp|stat", sc$stat); put("scan_lp|k", sc$k); put("scan_lp|m", sc$m)
sc <- rkhs_scan(x60, 0.25, kernel = "gaussian")
put("scan_rkhs|stat", sc$stat); put("scan_rkhs|k", sc$k)
put("scan_rkhs|m", sc$m)

for (g in c(0.25, 0.4)) {
  put(sprintf("far1|g=%s", g), far1_statistic(F4, g))
}
put("far1_beta3|g=0.4",
    far1_statistic(F4, 0.4, beta_t = far1_betas$beta3))
sc <- far1_scan(F4, 0.4)
put("scan_far1|stat", sc$stat); put("scan_far1|k", sc$k)
put("scan_far1|m", sc$m)

## ---- critical values (every shipped cache cell) -----------------------------
for (g in c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.495)) {
  put(sprintf("cv_wilcoxon|g=%s|b=0", g), wilcoxon_critical_value(g))
}
put("cv_wilcoxon|g=0.5|b=1", wilcoxon_critical_value(0.5, beta = 1))
put("cv_wilcoxon_s2|g=0.25", wilcoxon_critical_value(0.25, splits = 2))
for (g in c(0, 0.25, 0.45)) for (p in c(2L, 4L, 10L)) {
  put(sprintf("cv_lp|g=%s|p=%d", g, p), lp_critical_value(g, p))
}
rk <- list(c("gaussian", ""), c("min", "rnorm.sq"), c("min", "runif"),
           c("2max", "runif"), c("2max", "rbeta.1.9"))
for (g in c(0, 0.25, 0.45)) for (cell in rk) {
  put(sprintf("cv_rkhs|g=%s|k=%s|b=%s", g, cell[1], cell[2]),
      rkhs_critical_value(g, cell[1], cell[2]))
}
for (g in c(0, 0.25, 0.4, 0.45)) {
  put(sprintf("cv_far1|g=%s", g), far1_critical_value(g))
}

## ---- telco_far1 dataset conversion + end-to-end ------------------------------
for (nm in names(telco_far1)) {
  put(sprintf("telco|%s|sum", nm), sum(telco_far1[[nm]]))
}
put("telco|cl2|y[1,1]",   telco_far1$cl2[1, 1])
put("telco|cl2|y[37,28]", telco_far1$cl2[37, 28])
put("telco|cl9|y[50,14]", telco_far1$cl9[50, 14])
sc <- wilcoxon_scan(project_l2(t(telco_far1$cl2)), 0.25)
put("telco|scan_cl2|stat", sc$stat)
put("telco|scan_cl2|k", sc$k); put("telco|scan_cl2|m", sc$m)

## ---- projections and reconstruction on the F5 input -------------------------
write_vec(project_l2(F5), "expected_l2.csv")
write_mat(reconstruct(F5, "fs_interp", K = 4)$X, "expected_fsinterp.csv")
write_mat(reconstruct(F5, "fs_ls", K = 4)$X, "expected_fsls.csv")
write_mat(reconstruct(F5, "bspline_pen", K = 4, lambda = 1e-4)$X,
          "expected_bsplinepen.csv")

df <- do.call(rbind, rows)
write.csv(df, file.path(out_dir, "expected_scalars.csv"), row.names = FALSE)
cat(sprintf("wrote %d scalar expectations + inputs to %s\n", nrow(df), out_dir))
