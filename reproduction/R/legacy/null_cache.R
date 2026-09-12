#' Precomputed null distributions (vendored legacy loaders).
#'
#' The Monte-Carlo CSVs live in the repository's `data-raw/cache-src/`
#' directory (the caches rescued verbatim from the original pipeline).
#' Recomputing precise quantiles per (gamma, ...) cell takes weeks, so runs
#' read from these caches rather than calling their own MC routine.
#'
#'   Antisymmetric (Ch3): T2_gamma_beta.csv     (long: gamma, beta, x;
#'                                               |W(t)-W(s)| * rho(|t-s|(1-|t-s|)))
#'   Lp (Ch4):            dist_p2_ws_n200_t.csv (wide: dist.p2.ws.g<g>.n200.p<p>)
#'   RKHS (Ch5):          SD_h_n2000_dist_t.csv (wide: SD.<g>.h.<kernel>.n2000[.<baseline>])
#'   plus the corrected-family caches (*_corrected.csv) and finite nulls.
#'
#' Loaders cache per-CSV read in RAW_CACHE so repeated lookups are free.

RAW_CACHE <- new.env(parent = emptyenv())

# Captured at source-time: directory of THIS file. The previous version
# called `sys.frame(1)$ofile` from inside raw_root(), but by then the
# source() frame has unwound and the lookup falls back to getwd() -- which
# breaks when smoke.R is launched from the project root. Capturing once at
# load-time gives a stable anchor that survives mclapply forks.
.NULL_CACHE_DIR <- local({
  ofile <- NULL
  i <- sys.nframe() - 1L
  while (i >= 1L) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && nzchar(of)) { ofile <- of; break }
    i <- i - 1L
  }
  if (is.null(ofile)) {
    args <- commandArgs(trailingOnly = FALSE)
    hit <- sub("--file=", "", args[grepl("^--file=", args)])
    if (length(hit)) ofile <- hit[1]
  }
  if (!is.null(ofile)) dirname(normalizePath(ofile, mustWork = FALSE))
  else NA_character_
})

raw_root <- function() {
  env <- Sys.getenv("REPRO_CACHE_DIR", unset = NA)
  if (!is.na(env) && nzchar(env)) return(env)
  # This file lives at <repo>/reproduction/R/legacy/; the rescued caches at
  # <repo>/data-raw/cache-src/.
  cur <- if (!is.na(.NULL_CACHE_DIR)) .NULL_CACHE_DIR else getwd()
  for (i in 1:6) {
    cand <- file.path(cur, "data-raw", "cache-src")
    if (dir.exists(cand)) return(cand)
    cur <- dirname(cur)
  }
  stop("Could not locate data-raw/cache-src/; set REPRO_CACHE_DIR to override.")
}

raw_read <- function(filename) {
  if (exists(filename, envir = RAW_CACHE, inherits = FALSE)) {
    return(get(filename, envir = RAW_CACHE))
  }
  path <- file.path(raw_root(), filename)
  if (!file.exists(path)) {
    # Regenerated caches are written by their generators to
    # reproduction/results/ first and only hand-promoted to
    # data-raw/cache-src/ once sanity-checked; fall back so a regeneration
    # run works before the promotion.
    alt <- if (!is.na(.NULL_CACHE_DIR)) {
      file.path(dirname(dirname(.NULL_CACHE_DIR)), "results", filename)
    } else {
      file.path("results", filename)
    }
    if (file.exists(alt)) path <- alt else
      stop("Precomputed null not found: ", path, " (also tried ", alt, ")")
  }
  df <- read.csv(path, check.names = FALSE)
  assign(filename, df, envir = RAW_CACHE)
  df
}

#' Closed-form Kuiper law: the gamma = 0 Wilcoxon null.
#'
#' At gamma = 0 the weight is identically 1 and T_0 is exactly the range of a
#' standard Brownian bridge, whose distribution is Kuiper's V:
#'   P(V <= v) = 1 - sum_{k>=1} 2 (4 k^2 v^2 - 1) exp(-2 k^2 v^2).
#' Printed Table 3.1's gamma = 0 row matches this law at all nine quantiles
#' (slightly below it, as the grid discretisation of a sup predicts), while
#' the cached MC column is provably wrong: its 95% point (1.770) EXCEEDS the
#' continuous-limit value 1.7473, which a discretised max cannot do, and its
#' maximum over 30,000 draws (2.274) sits below the exact 99.9% point (2.303).
kuiper_cdf <- function(v) {
  vapply(v, function(vi) {
    if (vi <= 0.3) return(0)
    k <- 1:1000
    min(1, max(0, 1 - sum(2 * (4 * k^2 * vi^2 - 1) * exp(-2 * k^2 * vi^2))))
  }, numeric(1))
}

kuiper_quantile <- function(p) {
  vapply(p, function(pp) {
    pp <- min(max(pp, 1e-9), 1 - 1e-12)
    uniroot(function(v) kuiper_cdf(v) - pp, c(0.3, 6),
            tol = .Machine$double.eps^0.5)$root
  }, numeric(1))
}

#' Deterministic quantile grid standing in for a Monte-Carlo sample of the
#' Kuiper law: qK((i - 0.5)/N), i = 1..N. Empirical quantiles of this vector
#' reproduce the exact quantiles to O(1/N), so every existing consumer of
#' `load_wilcoxon_null` (which all take empirical quantiles) works unchanged.
#' Built by inverting a dense tabulated CDF (vectorised, ~0.1 s) rather than
#' 200k uniroot calls; interpolation error is ~1e-5, far below the 1e-3 the
#' 3-decimal critical values need. `kuiper_quantile` stays uniroot-exact for
#' spot checks against printed Table 3.1.
.kuiper_grid <- local({
  g <- NULL
  function(N = 200000L) {
    if (is.null(g)) {
      v   <- seq(0.3, 6, length.out = 20001L)
      vk  <- outer(v^2, (1:100)^2)          # v^2 * k^2 grid
      cdf <- 1 - rowSums(2 * (4 * vk - 1) * exp(-2 * vk))
      cdf <- cummax(pmin(1, pmax(0, cdf)))  # clamp + enforce monotonicity
      keep <- !duplicated(cdf)
      g <<- approx(cdf[keep], v[keep], xout = (seq_len(N) - 0.5) / N,
                   rule = 2)$y
    }
    g
  }
})

#' Antisymmetric-kernel (Wilcoxon) 2-sided limit distribution at (gamma, beta).
#' The CSV only varies beta at gamma = 0.5; elsewhere beta is forced to 0.
#' `exact_gamma0 = TRUE` serves the (0, 0) cell from the closed-form Kuiper
#' law instead of the cache (see `kuiper_cdf` above for why the cached column
#' is wrong); the gamma > 0 columns match printed Table 3.1 within MC noise
#' and stay cached either way.
load_wilcoxon_null <- function(gamma, beta = 0, exact_gamma0 = FALSE) {
  if (gamma != 0.5) beta <- 0
  if (isTRUE(exact_gamma0) && gamma == 0 && beta == 0) {
    return(.kuiper_grid())
  }
  df <- raw_read("T2_gamma_beta.csv")
  x <- df$x[df$gamma == gamma & df$beta == beta]
  if (!length(x)) {
    stop(sprintf("No Wilcoxon null for (gamma=%s, beta=%s). Available gammas: %s",
                 gamma, beta, paste(sort(unique(df$gamma)), collapse = ", ")))
  }
  x
}

#' Lp Wiener-sheet null at (gamma, p), n = 200.
load_lp_null <- function(gamma, p) {
  df <- raw_read("dist_p2_ws_n200_t.csv")
  col <- sprintf("dist.p2.ws.g%s.n200.p%d", format(gamma, trim = TRUE), p)
  if (!col %in% names(df)) {
    stop(sprintf("No Lp null column %s. Have: %s",
                 col, paste(names(df), collapse = ", ")))
  }
  df[[col]]
}

#' Birbilas's E_{q,gamma} asymptotic null. Indexed by gamma. Generated by
#' `generators/birbilas_null.R`, written to `.raw/E_q_gamma.csv`.
load_birbilas_null <- function(gamma) {
  df <- raw_read("E_q_gamma.csv")
  x <- df$x[df$gamma == gamma]
  if (!length(x)) {
    stop(sprintf("No Birbilas null for gamma=%s. Available: %s",
                 gamma, paste(sort(unique(df$gamma)), collapse = ", ")))
  }
  x
}

#' Birbilas's finite-sample (adjusted) null: empirical n^{-1/2 + gamma}
#' E_{n,beta,gamma} draws under H_0 at the target (n, gamma, beta_id).
#' Indexed by all three. Generated by `generators/birbilas_finite_null.R`,
#' written to `.raw/E_n_beta_gamma_finite.csv`.
#' `variant = "corrected"` reads the cache regenerated with
#' `birbilas_statistic_corrected` (printed length^gamma denominator); the two
#' statistics differ per window, so their finite nulls must never be mixed.
load_birbilas_finite_null <- function(n, gamma, beta_id = "beta1",
                                      variant = "legacy") {
  df <- raw_read(if (identical(variant, "corrected"))
                   "E_n_beta_gamma_finite_corrected.csv"
                 else "E_n_beta_gamma_finite.csv")
  x <- df$x[df$n == n & df$gamma == gamma & df$beta_id == beta_id]
  if (!length(x)) {
    keys <- unique(df[, c("n", "gamma", "beta_id")])
    stop(sprintf(
      "No Birbilas finite null for (n=%s, gamma=%s, beta_id=%s). Available cells: %s. Populate via generators/birbilas_finite_null.R.",
      n, gamma, beta_id,
      paste(apply(keys, 1L, function(r) sprintf("(n=%s,g=%s,%s)", r[1], r[2], r[3])),
            collapse = ", ")))
  }
  x
}

#' RKHS finite-sample null on F4 (FAR(1) + Wiener) data, indexed by
#' (n, gamma, projection, kernel). Generated by
#' `generators/rkhs_F4_null.R`. Kernel and projection are limited to a
#' hand-picked set of (projection, kernel) pairs (see the generator);
#' lookup errors list the available cells.
load_rkhs_F4_null <- function(n, gamma, projection, kernel_name,
                              variant = "legacy") {
  df <- raw_read(if (identical(variant, "corrected"))
                   "SD_h_F4_finite_corrected.csv"
                 else "SD_h_F4_finite.csv")
  x <- df$x[df$n == n & df$gamma == gamma &
            df$projection == projection & df$kernel == kernel_name]
  if (!length(x)) {
    keys <- unique(df[, c("n", "gamma", "projection", "kernel")])
    stop(sprintf(
      "No RKHS F4 null for (n=%s, gamma=%s, projection=%s, kernel=%s). Available cells: %s. Populate via generators/rkhs_F4_null.R.",
      n, gamma, projection, kernel_name,
      paste(apply(keys, 1L, function(r)
        sprintf("(n=%s,g=%s,%s,%s)", r[1], r[2], r[3], r[4])),
        collapse = ", ")))
  }
  x
}

#' RKHS null at (gamma, kernel, baseline), n = 2000.
#' `baseline` is one of "", "runif", "rnorm.sq" -- matching the CSV suffixes.
#' Only some (kernel, baseline) combinations were precomputed; see column names
#' in SD_h_n2000_dist_t.csv.
#' `variant = "corrected"` reads the cache regenerated with
#' `rkhs_statistic_corrected`; the legacy cache tabulates the biased statistic
#' (see methods/rkhs/run.R) and calibrates only that one.
load_rkhs_null <- function(gamma, kernel, baseline = "", variant = "legacy") {
  df <- raw_read(if (identical(variant, "corrected"))
                   "SD_h_n2000_corrected.csv"
                 else "SD_h_n2000_dist_t.csv")
  col <- sprintf("SD.%s.h.%s.n2000%s",
                 format(gamma, trim = TRUE),
                 kernel,
                 if (nzchar(baseline)) paste0(".", baseline) else "")
  if (!col %in% names(df)) {
    stop(sprintf("No RKHS null column %s. Have: %s",
                 col, paste(names(df), collapse = ", ")))
  }
  df[[col]]
}
