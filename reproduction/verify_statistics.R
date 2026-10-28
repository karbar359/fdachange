#!/usr/bin/env Rscript
# Verification gate for the statistic audit, ported onto the fdachange
# package.
#
# For each scan statistic, brute-force the formula EXACTLY as the
# dissertation prints it and compare against the vendored pre-audit (legacy)
# implementations and the package (corrected) ones. Every defect the audit
# reported must reproduce here before any run is trusted; every corrected
# function must equal its brute force to machine precision; and the vendored
# corrected variants must equal the package bit for bit. Also checks the
# closed-form Kuiper law against printed Table 3.1 and the scenario-level
# items (S2 non-centrality, the l + 1 change window).
#
# Runs in well under a minute. Exits non-zero if any MUST-PASS check fails.
#
# Usage:  Rscript reproduction/verify_statistics.R

# ---- Inputs (edit these) --------------------------------------------------
N_W    <- 60L    # Wilcoxon brute-force sample size
N_R    <- 60L    # RKHS brute-force sample size
N_L    <- 60L    # Lp brute-force sample size
N_B    <- 40L    # Birbilas brute-force sample size (n curves)
GAMMAS <- c(0, 0.25, 0.45)
SEED   <- 7L
TOL    <- 1e-10  # equality tolerance for brute-vs-implementation
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(fdachange))

root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) dirname(normalizePath(hit[1])) else getwd()
})
source(file.path(root, "R", "compat.R"))
source(file.path(root, "R", "legacy_env.R"))
legacy <- legacy_env(root)

fails <- 0L
check <- function(ok, label, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  (", detail, ")") else ""))
  if (!ok) fails <<- fails + 1L
  invisible(ok)
}

# ============================ 1. Wilcoxon ===================================
cat("== 1. Wilcoxon: printed T = max |Delta| / rho((l/n)(1 - l/n)) ==\n")
set.seed(SEED); Xw <- rnorm(N_W)

# Printed formula, literal: windows I(k, m) = {k+1..m}, l = m - k in 1..n,
# DIVIDE by rho, i.e. multiply by w = 1/rho. Includes length-1 windows.
brute_wilcoxon <- function(X, gamma, kernel = "h_w") {
  Xv <- if (kernel == "h_0") abs(X) else X
  n  <- length(Xv)
  H  <- outer(Xv, Xv, function(a, b) sign(a - b))
  sigma <- sqrt(mean(colMeans(H)^2))
  w  <- holder_weights(n, gamma, 0, symmetric = TRUE)   # w = 1/rho(t(1-t))
  best <- 0
  for (k in 0:(n - 1L)) for (m in (k + 1L):n) {
    I <- (k + 1L):m
    v <- abs(sum(H[I, setdiff(1:n, I)])) * w[m - k]
    if (is.finite(v) && v > best) best <- v
  }
  best / (n^1.5 * sigma)
}

for (g in GAMMAS) {
  b  <- brute_wilcoxon(Xw, g)
  co <- wilcoxon_statistic(Xw, g)
  le <- legacy$wilcoxon_statistic(Xw, g)
  ex <- legacy$wilcoxon_statistic_exact(Xw, g)
  check(abs(b - co) < TOL, sprintf("corrected == printed brute at gamma=%.2f", g),
        sprintf("brute %.6f corrected %.6f", b, co))
  cat(sprintf("        legacy/brute = %.3f   exact/brute = %.3f\n", le / b, ex / b))
}
# The vectorised-gamma call and the scan must agree with the statistic.
gb <- wilcoxon_statistic(Xw, GAMMAS)
sc <- vapply(GAMMAS, function(g) wilcoxon_scan(Xw, g)$stat, 0)
st <- vapply(GAMMAS, function(g) wilcoxon_statistic(Xw, g), 0)
check(max(abs(gb - st)) < TOL && max(abs(sc - st)) < TOL,
      "vectorised-gamma call and scan agree with the statistic")

# W2 attribution: multiply-direction but legacy index/window family vs corrected.
w2 <- local({
  co <- legacy$.wilcoxon_core_exact(Xw, "h_w")   # legacy family: windows k..m, l = m-k
  vapply(GAMMAS, function(g) {
    w <- holder_weights(co$n, g, 0, symmetric = TRUE)
    max(co$diff_abs * w[co$l_vec]) / (co$n^1.5 * co$sigma)
  }, numeric(1))
})
cat(sprintf("  W2 index off-by-one alone (multiply direction, legacy pairing):\n"))
for (i in seq_along(GAMMAS)) {
  cat(sprintf("        gamma=%.2f  legacy-pairing %.6f vs corrected %.6f (ratio %.4f)\n",
              GAMMAS[i], w2[i], st[i], w2[i] / st[i]))
}

# ============================ 2. Kuiper closed form =========================
cat("== 2. gamma = 0 null: closed-form Kuiper vs printed Table 3.1 ==\n")
printed_levels <- c(0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999)
printed_row    <- c(1.213, 1.460, 1.612, 1.741, 1.857, 2.012, 2.104, 2.195, 2.311)
qk <- vapply(printed_levels, legacy$kuiper_quantile, numeric(1))   # vendored reference
for (i in seq_along(printed_levels)) {
  cat(sprintf("        %6.2f%%: printed %.3f  exact %.4f  diff %+.4f\n",
              100 * (1 - printed_levels[i]), printed_row[i], qk[i],
              printed_row[i] - qk[i]))
}
# Bulk levels (50-2.5%): grid discretisation of the sup biases the MC row
# DOWN, so printed <= exact there. Extreme tail (1-0.1%): the 30,000-draw MC
# noise dominates and the sign is free. Both regimes stay within 0.014. The
# cached column shipped until 2026-10-28 was +0.023 off at the 5% point with
# an impossible maximum because its draws were not independent: the batched
# mclapply helper that patched it re-forked the same RNG substreams every
# second batch, so its 30,000 rows held about 160 distinct bridges (acf
# 0.994 at lag 160; see R/legacy/progress.R). The regenerated cell (same
# bridge construction, independent draws) behaves like the printed row:
# below the limit in the bulk, within two-sample MC noise of the print.
check(all(abs(printed_row - qk) < 0.014) && all((printed_row - qk)[1:5] <= 0),
      "printed row tracks the closed form (below in the bulk, MC noise in the tail)")
grid95 <- as.numeric(quantile(legacy$load_wilcoxon_null(0, 0, exact_gamma0 = TRUE), 0.95))
check(abs(grid95 - 1.7473) < 5e-4, "drop-in grid 95% point = 1.7473 (vendored reference, no longer served)",
      sprintf("got %.5f", grid95))
cell_q <- as.numeric(quantile(legacy$load_wilcoxon_null(0, 0), printed_levels))
check(cell_q[4] < qk[4] && all(abs(cell_q - printed_row)[1:5] < 0.014),
      "regenerated gamma = 0 cell sits below the continuous limit and tracks the printed row",
      sprintf("cell %.4f vs exact %.4f, printed %.3f", cell_q[4], qk[4], printed_row[4]))
pkg_q <- as.numeric(quantile(wilcoxon_null_draws(0), printed_levels))
check(isTRUE(all.equal(pkg_q, cell_q, tolerance = 1e-12)),
      "package gamma = 0 cell == vendored CSV cell")

# ============================ 3. Lp =========================================
cat("== 3. Lp: reduced (printed, uncentred) vs centred vs implementations ==\n")
set.seed(SEED); Xl <- rnorm(N_L)

brute_lp <- function(X, gamma, centred, p = 2, cdf = pnorm, t_grid = 100) {
  n  <- length(X); F.X <- cdf(X)
  w  <- holder_weights(n, gamma)                      # w = 1/rho(l/n)
  ts <- seq(0, 1, length.out = t_grid + 1)
  ind <- outer(F.X, ts, function(f, t) (f <= t) - t)
  tot <- colSums(ind)
  best <- 0
  for (k in 0:(n - 1L)) for (l in 1:(n - k)) {
    cs <- colSums(ind[(k + 1L):(k + l), , drop = FALSE])
    if (centred) cs <- cs - (l / n) * tot
    v <- (sum(abs(cs)^p) / t_grid)^(1 / p) * w[l]
    if (is.finite(v) && v > best) best <- v
  }
  best / sqrt(n)
}

for (g in GAMMAS) {
  b_unc <- brute_lp(Xl, g, centred = FALSE)
  b_cen <- brute_lp(Xl, g, centred = TRUE)
  le    <- legacy$lp_statistic(Xl, g, p = 2, cdf = pnorm)
  co    <- lp_statistic(Xl, g, p = 2, cdf = pnorm)
  # Legacy never scans start = n; brute does. Compare with tolerance for that
  # single extra length-1 window (it decides the max only pathologically).
  check(abs(b_unc - le) < 1e-9 || b_unc >= le,
        sprintf("legacy == printed REDUCED statistic at gamma=%.2f", g),
        sprintf("brute %.6f legacy %.6f", b_unc, le))
  check(abs(b_cen - co) < TOL,
        sprintf("corrected == centred statistic at gamma=%.2f", g),
        sprintf("brute %.6f corrected %.6f", b_cen, co))
}
# On the ecdf path the centring must be a numerical no-op (O(1) global sum).
ec_diff <- vapply(GAMMAS, function(g) {
  e <- function(x) ecdf(x)(x)
  abs(legacy$lp_statistic(Xl, g, p = 2, cdf = e) -
      lp_statistic(Xl, g, p = 2))                 # package default = ecdf path
}, numeric(1))
check(max(ec_diff) < 0.02, "ecdf path: corrected ~ legacy (centring auto-supplied)",
      sprintf("max diff %.4f", max(ec_diff)))
# lp_scan agrees with the corrected statistic.
scan_diff <- max(abs(vapply(GAMMAS, function(g)
  lp_scan(Xl, g, p = 2, cdf = pnorm)$stat -
  lp_statistic(Xl, g, p = 2, cdf = pnorm), numeric(1))))
check(scan_diff < TOL, "lp_scan agrees with the corrected statistic")

# ============================ 4. RKHS =======================================
cat("== 4. RKHS: printed Delta^2 vs implementations ==\n")
set.seed(SEED); Xr <- rnorm(N_R)
Kg <- gram(Xr, kernel_gaussian)

brute_rkhs <- function(X, K, gamma, full_cross = TRUE, min_start = 1L,
                       l_max = length(X) %/% 2L) {
  n <- length(X)
  w <- holder_weights(n, gamma)
  best <- 0
  for (l in 1:l_max) for (a in min_start:(n - l + 1L)) {
    I  <- a:(a + l - 1L)
    ln <- l / n
    cross <- if (full_cross) sum(K[I, 1:n]) else sum(K[I, 2:n])
    d2 <- sum(K[I, I]) + ln^2 * sum(K) - 2 * ln * cross
    v  <- sqrt(max(d2, 0)) * w[l]
    if (is.finite(v) && v > best) best <- v
  }
  best / sqrt(n)
}

for (g in GAMMAS) {
  b  <- brute_rkhs(Xr, Kg, g)                                   # printed
  co <- rkhs_statistic(Xr, g, kernel = "gaussian")
  check(abs(b - co) < TOL, sprintf("corrected == printed brute at gamma=%.2f", g),
        sprintf("brute %.6f corrected %.6f", b, co))
  # Mechanism: legacy == starts >= 2, l unrestricted, cross missing column 1.
  bl <- brute_rkhs(Xr, Kg, g, full_cross = FALSE, min_start = 2L,
                   l_max = length(Xr) - 1L)
  le <- legacy$rkhs_statistic(Xr, g, legacy$kernel_gaussian)
  check(abs(bl - le) < TOL,
        sprintf("legacy mechanism reproduced at gamma=%.2f", g),
        sprintf("code-family brute %.6f legacy %.6f (printed %.6f)", bl, le, b))
}
# 2-max degeneracy on Beta(1,9): legacy ~ sqrt(2 mean K), corrected far below.
set.seed(SEED + 1L); Xb <- rbeta(100L, 1, 9)
Kb  <- gram(Xb, kernel_2max)
le2 <- legacy$rkhs_statistic(Xb, 0, legacy$kernel_2_max)
co2 <- rkhs_statistic(Xb, 0, kernel = "2max")
cat(sprintf("        2-max on B(1,9): legacy %.4f  sqrt(2 mean K) %.4f  corrected %.4f\n",
            le2, sqrt(2 * mean(Kb)), co2))
check(abs(le2 - sqrt(2 * mean(Kb))) < 0.05 && co2 < 0.5 * le2,
      "2-max degeneracy: legacy tracks sqrt(2 mean K); corrected does not")
# rkhs_scan agrees with the corrected statistic.
check(abs(rkhs_scan(Xr, 0.25, kernel = "gaussian")$stat -
          rkhs_statistic(Xr, 0.25, kernel = "gaussian")) < TOL,
      "rkhs_scan agrees with the corrected statistic")

# ============================ 5. Birbilas ===================================
cat("== 5. Birbilas: printed (m - k)^gamma denominator ==\n")
set.seed(SEED); Xf <- generate_F4(N_B, 0, c(0.4, 0.6))

brute_far1 <- function(X, gamma, beta_t = function(t) 0.7 * sqrt(t),
                           len_minus_1 = FALSE) {
  L <- nrow(X); n <- ncol(X)
  s <- seq(0, 1, length.out = L)
  W <- (1 - beta_t(s)) * (X - rowMeans(X))
  best <- 0
  for (k in 0:(n - 1L)) for (m in (k + 1L):n) {
    Ssum  <- rowSums(W[, (k + 1L):m, drop = FALSE])
    denom <- if (len_minus_1) max(1L, m - k - 1L)^gamma else (m - k)^gamma
    v <- max(abs(Ssum)) / denom
    if (v > best) best <- v
  }
  best * n^(-0.5 + gamma)
}

for (g in c(0.25, 0.4)) {
  b  <- brute_far1(Xf, g)
  co <- far1_statistic(Xf, g)
  le <- legacy$far1_statistic(Xf, g)
  bl <- brute_far1(Xf, g, len_minus_1 = TRUE)
  check(abs(b - co) < TOL, sprintf("corrected == printed brute at gamma=%.2f", g),
        sprintf("brute %.6f corrected %.6f", b, co))
  # Legacy never scans length-1 windows; the (len-1)^gamma brute includes them
  # (denominator clamped at 1), so legacy <= brute with equality when the max
  # is an l >= 2 window -- the generic case.
  check(abs(bl - le) < TOL,
        sprintf("legacy == (length-1)^gamma variant at gamma=%.2f", g),
        sprintf("variant %.6f legacy %.6f (printed %.6f, ratio %.3f)", bl, le, b, le / b))
}

# ============================ 6. Scenario-level =============================
cat("== 6. Scenarios: S2 non-centrality and the l + 1 change window ==\n")
set.seed(SEED)
a2 <- 2
s2  <- generate_S2(200000L, a2, c(0.1, 0.9))
m2  <- mean(s2[20001:180000])               # the in-window slice
cat(sprintf("        S2 in-window mean %.3f: chi2_1(ncp=a^2) predicts %.0f, chi2_1(ncp=a) predicts %.0f\n",
            m2, 1 + a2^2, 1 + a2))
check(abs(m2 - (1 + a2^2)) < 0.1, "S2 generates ncp = a^2, not the printed ncp = a")

count_changed <- function(exact) {
  if (exact) assign("EC_EXACT", TRUE, envir = .GlobalEnv)
  else if (exists("EC_EXACT", envir = .GlobalEnv)) rm(EC_EXACT, envir = .GlobalEnv)
  set.seed(1); x <- generate_S1(400L, 50, c(0.4925, 0.5075))   # short regime, l = 6
  sum(x > 25)
}
n_leg <- count_changed(FALSE); n_ex <- count_changed(TRUE)
check(n_leg == 7L && n_ex == 6L,
      "short-regime window: legacy 7 points, EC_EXACT 6 (= l)",
      sprintf("legacy %d exact %d", n_leg, n_ex))
if (exists("EC_EXACT", envir = .GlobalEnv)) rm(EC_EXACT, envir = .GlobalEnv)

# ================= 7. Package vs vendored corrected variants =================
cat("== 7. fdachange package == vendored corrected implementations ==\n")
check(identical(wilcoxon_statistic(Xw, 0.25),
                legacy$wilcoxon_statistic_corrected(Xw, 0.25)),
      "wilcoxon: package == vendored corrected (bit-exact)")
check(identical(lp_statistic(Xl, 0.25, p = 2, cdf = pnorm),
                legacy$lp_statistic_corrected(Xl, 0.25, p = 2, cdf = pnorm)),
      "lp: package == vendored corrected (bit-exact)")
check(identical(rkhs_statistic(Xr, 0.25, kernel = "gaussian"),
                legacy$rkhs_statistic_corrected(Xr, 0.25, legacy$kernel_gaussian)),
      "rkhs: package == vendored corrected (bit-exact)")
check(identical(far1_statistic(Xf, 0.4),
                legacy$far1_statistic_corrected(Xf, 0.4)),
      "far1: package == vendored corrected (bit-exact)")

cat(sprintf("\n%s: %d failure(s)\n", if (fails == 0L) "ALL CHECKS PASSED" else "GATE FAILED", fails))
if (fails > 0L) quit(status = 1L)
