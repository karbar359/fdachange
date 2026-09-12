#!/usr/bin/env Rscript
# Monte Carlo calibration demo: the thesis's per-(n, gamma) remedy, end to end.
#
# Demonstrates, with the pipeline's own statistic code, that
#   (1) PROBLEM: the asymptotic critical values are conservative at small n
#       for the Wilcoxon scan -- size far below alpha and power collapse at
#       gamma = 0.45 -- and anticonservative for the Lp scan at gamma = 0
#       (over-rejection);
#   (2) FIX: Monte Carlo calibration at the actual (n, gamma) -- simulate the
#       null at sample size n, take the empirical (1-alpha) quantile as the
#       critical value -- restores exact size and the lost power.
#
# The calibrated table is still DATA-FREE: both statistics are rank
# statistics under continuous F (Wilcoxon via ranks; Lp after the PIT), so
# one H0 simulation per (n, gamma) -- run under any continuous law -- is
# valid for every continuous data distribution. This is the "per (n, gamma)
# but data-free" remedy named in the thesis conclusions.
#
# Design per n:
#   batch A (N_CAL H0 samples)  -> calibrated CV = order statistic at rank
#                                  ceil((N_CAL+1)(1-alpha)) per gamma
#   batch B (N_REP H0 samples)  -> size under asymptotic CV vs calibrated CV
#   batch C (N_REP H1 samples)  -> power under both CVs (Wilcoxon only;
#                                  S1-style mean shift +A on the centered
#                                  window of length n*SEG_FRAC, the thesis's
#                                  long-segment design)
# Batches are independent; all gammas share one statistic evaluation per
# sample (batched implementations).
#
# Usage: Rscript reproduction/mc_calibration_demo.R
# REPRO_CACHE_DIR may be set to override the cache lookup (asymptotic
# critical values come from the same cached nulls the thesis pipeline uses).

# ---- Inputs (edit these) --------------------------------------------------
if (!exists("QUICK")) QUICK <- TRUE   # TRUE: sub-minute sanity run; FALSE: full
N_S      <- if (QUICK) c(100L, 200L) else c(100L, 200L, 400L)
N_CAL    <- if (QUICK) 300L  else 2000L   # H0 draws for calibration (batch A)
N_REP    <- if (QUICK) 200L  else 1000L   # evaluation replicates (batches B, C)
GAMMA_S  <- c(0, 0.25, 0.45)
ALPHA    <- 0.05
A_SHIFT  <- 1                      # S1-style mean shift on the changed window
SEG_FRAC <- 1 / 5                  # window length as fraction of n (centered)
LP_TGRID <- 100L                   # t-grid of the Lp statistic (pipeline value)
NCORES   <- max(1L, parallel::detectCores() - 2L)
SEED     <- 42L
OUT_CSV  <- NULL                   # NULL -> data/mc_calibration_demo.csv on a
                                   # full run; QUICK runs write to results/
                                   # (gitignored) so they never clobber the
                                   # committed full-run table.
# ---------------------------------------------------------------------------

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- sub("--file=", "", args[grepl("^--file=", args)])
  if (length(hit)) return(normalizePath(hit[1]))
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) return(normalizePath(ofile))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(script_path)
stopifnot(dir.exists(file.path(root, "R", "legacy")))

source(file.path(root, "R", "legacy", "integral_image.R"))
source(file.path(root, "R", "legacy", "weights.R"))
source(file.path(root, "R", "legacy", "null_cache.R"))
source(file.path(root, "R", "legacy", "progress.R"))
source(file.path(root, "R", "legacy", "method_wilcoxon.R"))
source(file.path(root, "R", "legacy", "method_lp.R"))

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

if (is.null(OUT_CSV)) {
  OUT_CSV <- file.path(root, if (QUICK) "results" else "data",
                       "mc_calibration_demo.csv")
}

# Lp statistic for all gammas in one pass over the (k, l) window loop.
# Mirrors methods/lp/run.R::lp_statistic exactly (same loops, same window
# convention); only the gamma-specific weight multiplication is batched.
# Cross-checked below against lp_statistic before use.
lp_statistic_gammas_local <- function(X, gammas, p = 2,
                                      cdf = function(x) pnorm(x, 0, 1),
                                      t_grid = 100) {
  n <- length(X)
  F.X <- cdf(X)
  W <- vapply(gammas, function(g) weight_cache(n, g), numeric(n))
  ts <- seq(0, 1, length.out = t_grid + 1)
  ind <- outer(F.X, ts, FUN = function(f, t) (f <= t) - t)
  best <- rep(0, length(gammas))
  for (k in seq_len(n - 1)) {
    cs <- ind[k, ]
    for (l in seq_len(n - k)) {
      if (l > 1) cs <- cs + ind[k + l - 1, ]
      s <- (sum(abs(cs)^p) / t_grid)^(1 / p)
      v <- s * W[l, ]
      upd <- v > best
      if (any(upd)) best[upd] <- v[upd]
    }
  }
  best / sqrt(n)
}

cat("# cross-check: batched Lp variant vs methods/lp/run.R::lp_statistic\n")
xchk <- rnorm(40)
ref  <- vapply(GAMMA_S, function(g) lp_statistic(xchk, g, p = 2, t_grid = LP_TGRID),
               numeric(1))
got  <- lp_statistic_gammas_local(xchk, GAMMA_S, p = 2, t_grid = LP_TGRID)
stopifnot(isTRUE(all.equal(ref, got)))
cat("#   ok (max abs diff = ", format(max(abs(ref - got)), digits = 3), ")\n", sep = "")

# One evaluation batch: nrep samples of size n, statistic for every gamma.
# shift = 0 -> H0; shift > 0 -> mean shift on the centered SEG_FRAC window.
run_batch <- function(nrep, n, shift, stat_fn, label) {
  l   <- as.integer(round(n * SEG_FRAC))
  k0  <- as.integer(floor(n / 2) - floor(l / 2))
  idx <- (k0 + 1L):(k0 + l)
  cat(sprintf("  %s (nrep=%d)\n", label, nrep))
  rows <- mclapply_with_progress(seq_len(nrep), function(i) {
    X <- rnorm(n)
    if (shift != 0) X[idx] <- X[idx] + shift
    stat_fn(X)
  }, ncores = NCORES)
  do.call(rbind, rows)
}

# Calibrated CV: order statistic at rank ceil((N+1)(1-alpha)) of the H0
# draws -- the exact-level permutation-test convention.
cal_cv <- function(stats_mat) {
  r <- ceiling((nrow(stats_mat) + 1L) * (1 - ALPHA))
  stopifnot(r <= nrow(stats_mat))
  apply(stats_mat, 2L, function(s) sort(s)[r])
}

# Corrected scan (2026-09 audit): the legacy statistic applied the weight in
# the reciprocal direction, which is what produced this demo's original
# dramatic asymptotic-vs-calibrated gap (cv 3.85 vs finite-n 0.85 at
# gamma = 0.45) -- most of that gap was the inversion, not finite-n slowness.
wilc_fn <- function(X) wilcoxon_statistic_gammas_corrected(X, GAMMA_S)
lp_fn   <- function(X) lp_statistic_gammas_local(X, GAMMA_S, p = 2, t_grid = LP_TGRID)

results <- list()
t_start <- Sys.time()

for (n in N_S) {
  # ---- Wilcoxon: problem + fix -------------------------------------------
  cat(sprintf("== wilcoxon, n = %d ==\n", n))
  cv_asym <- vapply(GAMMA_S, function(g)
                      wilcoxon_critical_value(n, g, alpha = ALPHA,
                                              exact_gamma0 = TRUE),
                    numeric(1))
  A <- run_batch(N_CAL, n, 0,       wilc_fn, "batch A: calibration (H0)")
  B <- run_batch(N_REP, n, 0,       wilc_fn, "batch B: size (H0)")
  C <- run_batch(N_REP, n, A_SHIFT, wilc_fn, "batch C: power (H1)")
  cv_mc <- cal_cv(A)
  for (j in seq_along(GAMMA_S)) {
    results[[length(results) + 1L]] <- data.frame(
      test = "wilcoxon", n = n, gamma = GAMMA_S[j],
      cv_asym = cv_asym[j], cv_mc = cv_mc[j],
      size_asym  = mean(B[, j] > cv_asym[j]), size_mc  = mean(B[, j] > cv_mc[j]),
      power_asym = mean(C[, j] > cv_asym[j]), power_mc = mean(C[, j] > cv_mc[j]))
  }

  # ---- Lp: the anticonservative mirror (size only) ------------------------
  cat(sprintf("== lp (p = 2), n = %d ==\n", n))
  cv_asym <- vapply(GAMMA_S, function(g)
    critical_value(lp_null_distribution(n, g, p = 2), ALPHA, 1), numeric(1))
  A <- run_batch(N_CAL, n, 0, lp_fn, "batch A: calibration (H0)")
  B <- run_batch(N_REP, n, 0, lp_fn, "batch B: size (H0)")
  cv_mc <- cal_cv(A)
  for (j in seq_along(GAMMA_S)) {
    results[[length(results) + 1L]] <- data.frame(
      test = "lp", n = n, gamma = GAMMA_S[j],
      cv_asym = cv_asym[j], cv_mc = cv_mc[j],
      size_asym  = mean(B[, j] > cv_asym[j]), size_mc  = mean(B[, j] > cv_mc[j]),
      power_asym = NA_real_, power_mc = NA_real_)
  }
}

res <- do.call(rbind, results)
res$alpha <- ALPHA; res$a <- A_SHIFT; res$seg_frac <- SEG_FRAC
res$n_cal <- N_CAL; res$n_rep <- N_REP; res$seed <- SEED; res$quick <- QUICK

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(res, OUT_CSV, row.names = FALSE)

cat(sprintf("\n# elapsed: %s   (QUICK = %s, NCORES = %d)\n",
            fmt_duration(as.numeric(difftime(Sys.time(), t_start, units = "secs"))),
            QUICK, NCORES))
cat(sprintf("# written: %s\n\n", OUT_CSV))
cat(sprintf("%-9s %5s %6s %8s %7s %10s %8s %11s %9s\n",
            "test", "n", "gamma", "cv_asym", "cv_mc",
            "size_asym", "size_mc", "power_asym", "power_mc"))
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-9s %5d %6.2f %8.3f %7.3f %10.3f %8.3f %11s %9s\n",
              res$test[i], res$n[i], res$gamma[i],
              res$cv_asym[i], res$cv_mc[i],
              res$size_asym[i], res$size_mc[i],
              ifelse(is.na(res$power_asym[i]), "-", sprintf("%.3f", res$power_asym[i])),
              ifelse(is.na(res$power_mc[i]),  "-", sprintf("%.3f", res$power_mc[i]))))
}
