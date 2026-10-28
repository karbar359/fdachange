#!/usr/bin/env Rscript
# One-time capture of golden values from the ORIGINAL thesis codebase
# (corrected family), written at full double precision (%.17g) to
# tests/testthat/testdata/golden_values.csv. The package test suite
# recomputes every value with package functions; agreement at 1e-12 is the
# permanent tie between the package and the audited implementation.
#
# Needs the thesis repo (PHD_NOTES_DIR env var or the default path below)
# with EnvStats installed, and the null caches in data-raw/cache-src.
# Run from the package root:  Rscript data-raw/capture_golden_values.R

NOTES <- Sys.getenv("PHD_NOTES_DIR")
if (!nzchar(NOTES)) {
  stop("set PHD_NOTES_DIR to the root of the original thesis repo")
}
OLD_CODE <- file.path(NOTES, "code")
Sys.setenv(PHD_RAW_DIR = normalizePath(file.path("data-raw", "cache-src")))

old <- new.env(parent = globalenv())
suppressPackageStartupMessages({
  library(EnvStats)   # old scenario files need it at source time
  for (f in c("shared/integral_image.R", "shared/weights.R",
              "shared/kernels.R", "shared/null_cache.R",
              "shared/projections.R",
              "methods/wilcoxon/run.R", "methods/lp/run.R",
              "methods/rkhs/run.R", "methods/birbilas/run.R",
              "scenarios/S1/data.R", "scenarios/S2/data.R",
              "scenarios/S3/data.R", "scenarios/S4/data.R",
              "scenarios/S5/data.R", "scenarios/F1/data.R",
              "scenarios/F2/data.R", "scenarios/F3/data.R",
              "scenarios/F4/data.R")) {
    sys.source(file.path(OLD_CODE, f), envir = old, keep.source = FALSE)
  }
})

rows <- list()
put <- function(id, value) {
  rows[[length(rows) + 1L]] <<- data.frame(id = id,
                                           value = sprintf("%.17g", value))
}

SEED <- 7L
GAMMAS <- c(0, 0.25, 0.45)
NS <- c(40L, 60L, 120L)

## ---- statistics ------------------------------------------------------------
for (n in NS) {
  set.seed(SEED + n); x <- rnorm(n)
  for (g in GAMMAS) {
    for (kern in c("h_w", "h_0")) {
      put(sprintf("wilcoxon|n=%d|g=%s|k=%s", n, g, kern),
          old$wilcoxon_statistic_corrected(x, g, kernel = kern))
    }
    for (p in c(2L, 4L)) {
      put(sprintf("lp|n=%d|g=%s|p=%d|cdf=pnorm", n, g, p),
          old$lp_statistic_corrected(x, g, p = p, cdf = pnorm))
      put(sprintf("lp|n=%d|g=%s|p=%d|cdf=ecdf", n, g, p),
          old$lp_statistic_corrected(x, g, p = p,
                                     cdf = function(z) ecdf(z)(z)))
    }
    for (kern in c("gaussian", "min", "2_max")) {
      put(sprintf("rkhs|n=%d|g=%s|k=%s", n, g, sub("_", "", kern)),
          old$rkhs_statistic_corrected(x, g, old$kernel_fn(kern)))
    }
  }
}
for (n in c(40L, 60L)) {
  set.seed(SEED + n); Xf <- old$generate_F4(n, 0, c(0.4, 0.6))
  for (g in c(0.25, 0.4)) {
    for (b in c("beta1", "beta3")) {
      put(sprintf("far1|n=%d|g=%s|b=%s|a=0", n, g, b),
          old$birbilas_statistic_corrected(Xf, g,
            beta_t = old$birbilas_beta_registry[[b]]))
    }
  }
  set.seed(SEED + n); Xf2 <- old$generate_F4(n, 2, c(0.25, 0.45))
  put(sprintf("far1|n=%d|g=0.4|b=beta1|a=2", n),
      old$birbilas_statistic_corrected(Xf2, 0.4))
}

## ---- scans (stat, k, m) ----------------------------------------------------
set.seed(SEED + 60L); x60 <- rnorm(60)
sc <- old$wilcoxon_scan_corrected(x60, 0.25)
put("scan|wilcoxon|stat", sc$stat); put("scan|wilcoxon|k", sc$k)
put("scan|wilcoxon|m", sc$m)
sc <- old$lp_scan(x60, 0.25, cdf = pnorm)
put("scan|lp|stat", sc$stat); put("scan|lp|k", sc$k); put("scan|lp|m", sc$m)
sc <- old$rkhs_scan(x60, 0.25, old$kernel_gaussian)
put("scan|rkhs|stat", sc$stat); put("scan|rkhs|k", sc$k)
put("scan|rkhs|m", sc$m)

## ---- critical values from every shipped cache cell -------------------------
wc <- read.csv(file.path("data-raw", "cache-src", "T2_gamma_beta.csv"))
for (i in seq_len(nrow(unique(wc[, c("gamma", "beta")])))) {
  cell <- unique(wc[, c("gamma", "beta")])[i, ]
  g <- cell$gamma; b <- cell$beta
  draws <- old$load_wilcoxon_null(g, b)
  for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
    put(sprintf("cv|wilcoxon|g=%s|b=%s|a=%s|s=%d", g, b, al, m),
        old$critical_value(draws, al, m))
  }
}
for (g in GAMMAS) for (p in c(2L, 4L, 10L)) {
  for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
    put(sprintf("cv|lp|g=%s|p=%d|a=%s|s=%d", g, p, al, m),
        old$critical_value(old$load_lp_null(g, p), al, m))
  }
}
rk_cells <- list(c("gaussian", ""), c("min", "rnorm.sq"), c("min", "runif"),
                 c("2.max", "runif"), c("2.max", "rbeta.1.9"))
for (g in GAMMAS) for (cell in rk_cells) {
  draws <- old$load_rkhs_null(g, cell[1], cell[2], variant = "corrected")
  for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
    put(sprintf("cv|rkhs|g=%s|k=%s|b=%s|a=%s|s=%d",
                g, sub(".", "", cell[1], fixed = TRUE), cell[2], al, m),
        old$critical_value(draws, al, m))
  }
}
for (g in c(0, 0.25, 0.4, 0.45)) {
  for (al in c(0.05, 0.01)) for (m in c(1L, 2L)) {
    put(sprintf("cv|far1|g=%s|a=%s|s=%d", g, al, m),
        old$critical_value(old$load_birbilas_null(g), al, m))
  }
}

## ---- scenario draw streams (sum pins the whole RNG stream) -----------------
assign("EC_EXACT", TRUE, envir = old)
for (s in c("S1", "S2", "S3", "S4", "S5")) {
  gen <- old[[paste0("generate_", s)]]
  set.seed(11L); put(sprintf("stream|%s|a=1|exact", tolower(s)),
                     sum(gen(50, a = 1)))
}
for (s in c("F1", "F3")) {
  gen <- old[[paste0("generate_", s)]]
  set.seed(11L); put(sprintf("stream|%s|a=3|exact", tolower(s)),
                     sum(gen(20, a = 3, l = 20)))
}
set.seed(11L); put("stream|f2|a=0.5|exact", sum(old$generate_F2(20, a = 0.5, l = 20)))
set.seed(11L); put("stream|f4|a=2", sum(old$generate_F4(20, a = 2, l = 20)))
rm("EC_EXACT", envir = old)
for (s in c("S1", "S3")) {
  gen <- old[[paste0("generate_", s)]]
  set.seed(11L); put(sprintf("stream|%s|a=1|legacy", tolower(s)),
                     sum(gen(50, a = 1)))
}

## ---- projections ------------------------------------------------------------
set.seed(11L); Xp <- old$generate_F1(20, a = 0, l = 20)
put("proj|l2|sum", sum(old$project_l2(Xp)))
set.seed(13L); put("proj|ito|sum", sum(old$project_ito(Xp)))
set.seed(13L); put("proj|c_max_var|sum", sum(old$project_c_max_var(Xp)))
set.seed(13L); put("proj|wiener|sum", sum(old$project_wiener(Xp)))
set.seed(13L); put("proj|pc1|sumabs", sum(abs(old$project_pc1(Xp))))

df <- do.call(rbind, rows)
out <- file.path("tests", "testthat", "testdata")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(df, file.path(out, "golden_values.csv"), row.names = FALSE)
cat(sprintf("wrote %s (%d values)\n",
            file.path(out, "golden_values.csv"), nrow(df)))
