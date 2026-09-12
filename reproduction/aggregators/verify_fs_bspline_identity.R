#!/usr/bin/env Rscript
# Verify the claim that bridges the thesis text and the implementation:
#
#   span{ order-2 B-splines, breaks k/2^K }  ==  span{ Faber-Schauder truncated
#                                                      at level K }
#
# and separate that (true) span identity from the (false) claim that the two
# induce the same MAP into that space.
#
# The thesis describes Faber-Schauder truncation
# (chapters/ch1_introduction/subsec_faber_schauder.tex); the code builds an
# order-2 B-spline basis (shared/projections.R:26, shared/reconstruction.R).
# The thesis rationale states "that is WHY FS and not e.g. B-splines"
# (stated in the thesis defence materials), which is only consistent with the
# code because of the identity checked here.
#
# Usage:  Rscript reproduction/verify_fs_bspline_identity.R

# ---- Inputs (edit these) --------------------------------------------------
K_GRID   <- 2:6      # dyadic truncation levels to check
N_SAMPLE <- 168L     # sampling-grid size used by the application
TOL      <- 1e-10    # residual tolerance for "same span"
SEED     <- 42L
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(fda) })

root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) dirname(dirname(normalizePath(hit[1]))) else getwd()
})
source(file.path(root, "R", "legacy", "reconstruction.R"))

set.seed(SEED)

#' Faber-Schauder system on [0,1] truncated at level K, exactly as defined in
#' subsec_faber_schauder.tex: the constant, the identity, and the dyadic hat
#' functions at levels j = 0..K-1. Dimension 2 + (2^K - 1) = 2^K + 1.
fs_system <- function(t, K) {
  cols <- list(rep(1, length(t)), t)
  for (j in 0:(K - 1L)) {
    for (k in 0:(2^j - 1L)) {
      left <- k / 2^j; mid <- (k + 0.5) / 2^j; right <- (k + 1) / 2^j
      hat <- pmax(0, pmin((t - left) / (mid - left), (right - t) / (right - mid)))
      cols[[length(cols) + 1L]] <- hat
    }
  }
  do.call(cbind, cols)
}

cat("== 1. Dimension and span identity ==\n")
ok_all <- TRUE
for (K in K_GRID) {
  t  <- seq(0, 1, length.out = N_SAMPLE)
  B  <- eval.basis(t, fs_basis(K))   # order-2 B-splines, dyadic breaks
  Fm <- fs_system(t, K)              # thesis Faber-Schauder system

  dim_ok <- (ncol(B) == 2^K + 1L) && (ncol(Fm) == 2^K + 1L)

  # Each basis must reproduce the other exactly (least squares residual ~ 0).
  r_bf <- max(abs(Fm - B  %*% qr.solve(B,  Fm)))   # FS in terms of B-splines
  r_fb <- max(abs(B  - Fm %*% qr.solve(Fm, B )))   # B-splines in terms of FS
  ok   <- dim_ok && r_bf < TOL && r_fb < TOL
  ok_all <- ok_all && ok

  cat(sprintf("  K=%d  dim(Bspline)=%2d  dim(FS)=%2d  (2^K+1=%2d)  resid FS<-B %.2e  B<-FS %.2e  %s\n",
              K, ncol(B), ncol(Fm), 2^K + 1L, r_bf, r_fb,
              if (ok) "OK" else "FAIL"))
}

cat("\n== 2. Same space, different maps ==\n")
# A rough curve, in the spirit of hourly balancing prices.
t   <- seq(0, 1, length.out = N_SAMPLE)
x   <- sin(6 * pi * t) + 0.4 * rnorm(N_SAMPLE)
X   <- matrix(x, ncol = 1L)
for (K in K_GRID) {
  a <- reconstruct(X, "fs_interp", K, grid = "input")$X
  b <- reconstruct(X, "fs_ls",     K, grid = "input")$X
  cat(sprintf("  K=%d  max|fs_interp - fs_ls| = %7.4f   (sd of curve = %.4f, ratio %.2f)\n",
              K, max(abs(a - b)), sd(x), max(abs(a - b)) / sd(x)))
}
cat("  Both land in the same span; they are different projections onto it.\n")
cat("  fs_interp reproduces the thesis formula, fs_ls is what Data2fd does.\n")

cat("\n== 3. FS truncation really is polygonal interpolation ==\n")
# The classical Faber-Schauder coefficients of a function f are
#   c0 = f(0),  c1 = f(1) - f(0),
#   lambda_{j,k} = f(mid) - (f(left) + f(right)) / 2,
# and the series truncated at level K equals the piecewise-linear interpolant
# of f through the 2^K + 1 dyadic points. Checked on an explicit f so that f
# can be evaluated exactly at the dyadic points (a sampled curve would force an
# interpolation step and beg the question).
f_true <- function(u) sin(6 * pi * u) + 0.3 * cos(11 * pi * u)
fine   <- seq(0, 1, length.out = 4097L)
for (K in K_GRID) {
  Fm <- fs_system(fine, K)
  coefs <- numeric(ncol(Fm))
  coefs[1] <- f_true(0)
  coefs[2] <- f_true(1) - f_true(0)
  idx <- 2L
  for (j in 0:(K - 1L)) {
    for (k in 0:(2^j - 1L)) {
      left <- k / 2^j; mid <- (k + 0.5) / 2^j; right <- (k + 1) / 2^j
      idx <- idx + 1L
      coefs[idx] <- f_true(mid) - (f_true(left) + f_true(right)) / 2
    }
  }
  series <- as.numeric(Fm %*% coefs)
  dy     <- seq(0, 1, by = 2^(-as.integer(K)))
  interp <- approx(dy, f_true(dy), xout = fine)$y
  cat(sprintf("  K=%d  max|FS truncation - polygonal interpolation| = %.3e\n",
              K, max(abs(series - interp))))
}

cat("\n== 4. The figure-code basis is one dimension short ==\n")
# archive/code/figures/regenerate_fsb.R:20 and code/plot_eysm_litgrid_motivation.R:90
# seed the list with the constant only and loop over ODD k, omitting psi(t) = t.
figure_basis <- function(t, N) {
  b <- list(rep(1, length(t)))
  for (n in 1:N) {
    s <- 1 / 2^n
    for (kk in seq(1, 2^n - 1, by = 2)) {
      c0 <- kk / 2^n
      b[[length(b) + 1L]] <- pmax(0, 1 - abs(t - c0) / s)
    }
  }
  do.call(cbind, b)
}
for (K in K_GRID) {
  fb <- figure_basis(t, K)
  cat(sprintf("  K=%d  figure-code dim = %2d, correct dim = %2d  (%s)\n",
              K, ncol(fb), 2^K + 1L,
              if (ncol(fb) == 2^K + 1L) "OK" else "MISSING the psi(t)=t term"))
}

cat(sprintf("\nSpan identity: %s\n", if (ok_all) "VERIFIED" else "FAILED"))
if (!ok_all) quit(status = 1L)
