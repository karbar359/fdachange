#!/usr/bin/env Rscript
# Convert the rescued null-cache CSVs (data-raw/cache-src/, copied from the
# thesis repo) into the compressed extdata lists the package reads. Asserts
# element-wise identical() between each CSV cell and the rds entry before
# writing, so critical values stay bit-identical to the original pipeline.
#
# Run from the package root:  Rscript data-raw/build_null_caches.R

src <- file.path("data-raw", "cache-src")
out <- file.path("inst", "extdata")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

fmt <- function(x) format(x, trim = TRUE, scientific = FALSE)

stopifnot(dir.exists(src))

## ---- Wilcoxon: T2_gamma_beta.csv (long: x, gamma, beta) -------------------
t2 <- read.csv(file.path(src, "T2_gamma_beta.csv"))
cells <- unique(t2[, c("gamma", "beta")])
wilcoxon <- list()
for (i in seq_len(nrow(cells))) {
  g <- cells$gamma[i]; b <- cells$beta[i]
  key <- sprintf("g%s_b%s", fmt(g), fmt(b))
  wilcoxon[[key]] <- t2$x[t2$gamma == g & t2$beta == b]
  stopifnot(identical(wilcoxon[[key]], t2$x[t2$gamma == g & t2$beta == b]))
}
cat(sprintf("wilcoxon: %d cells (%s)\n", length(wilcoxon),
            paste(names(wilcoxon), collapse = ", ")))
saveRDS(wilcoxon, file.path(out, "wilcoxon_null.rds"), compress = "xz")

## ---- Lp: dist_p2_ws_n200_t.csv (wide: dist.p2.ws.g<g>.n200.p<p>) ----------
lpdf <- read.csv(file.path(src, "dist_p2_ws_n200_t.csv"), check.names = FALSE)
lp <- list()
for (col in names(lpdf)) {
  m <- regmatches(col, regexec("^dist\\.p2\\.ws\\.g([0-9.]+)\\.n200\\.p([0-9]+)$", col))[[1]]
  stopifnot(length(m) == 3)
  key <- sprintf("g%s_p%s", fmt(as.numeric(m[2])), m[3])
  lp[[key]] <- lpdf[[col]]
  stopifnot(identical(lp[[key]], lpdf[[col]]))
}
cat(sprintf("lp: %d cells (%s)\n", length(lp), paste(names(lp), collapse = ", ")))
saveRDS(lp, file.path(out, "lp_null.rds"), compress = "xz")

## ---- RKHS: SD_h_n2000_corrected.csv (wide: SD.<g>.h.<kernel>.n2000[.<b>]) -
rkdf <- read.csv(file.path(src, "SD_h_n2000_corrected.csv"), check.names = FALSE)
rkhs <- list()
for (col in names(rkdf)) {
  m <- regmatches(col, regexec(
    "^SD\\.([0-9.]+)\\.h\\.(gaussian|min|2\\.max)\\.n2000(\\.(.+))?$", col))[[1]]
  stopifnot(length(m) >= 3)
  kernel <- sub("2.max", "2max", m[3], fixed = TRUE)
  baseline <- if (length(m) >= 5 && nzchar(m[5])) m[5] else ""
  key <- sprintf("g%s_%s%s", fmt(as.numeric(m[2])), kernel,
                 if (nzchar(baseline)) paste0("_", baseline) else "")
  rkhs[[key]] <- rkdf[[col]]
  stopifnot(identical(rkhs[[key]], rkdf[[col]]))
}
cat(sprintf("rkhs: %d cells (%s)\n", length(rkhs), paste(names(rkhs), collapse = ", ")))
saveRDS(rkhs, file.path(out, "rkhs_null.rds"), compress = "xz")

## ---- Birbilas: E_q_gamma.csv (long: x, gamma) -----------------------------
bbdf <- read.csv(file.path(src, "E_q_gamma.csv"))
birbilas <- list()
for (g in unique(bbdf$gamma)) {
  key <- sprintf("g%s", fmt(g))
  birbilas[[key]] <- bbdf$x[bbdf$gamma == g]
  stopifnot(identical(birbilas[[key]], bbdf$x[bbdf$gamma == g]))
}
cat(sprintf("birbilas: %d cells (%s)\n", length(birbilas),
            paste(names(birbilas), collapse = ", ")))
saveRDS(birbilas, file.path(out, "birbilas_null.rds"), compress = "xz")

cat("\nextdata sizes:\n")
for (f in list.files(out, full.names = TRUE)) {
  cat(sprintf("  %-24s %8.0f KB\n", basename(f), file.size(f) / 1024))
}
