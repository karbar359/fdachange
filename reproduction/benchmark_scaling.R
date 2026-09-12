#!/usr/bin/env Rscript
# Minimal limit test: time and memory of the three scan statistics vs n.
#
# Motivated by a panel-report question on the tests' time/memory
# complexity and optimization headroom). Measures, per sample size n:
#   scan       : the shipped statistics exactly as the simulations ran them
#                (methods/{wilcoxon,rkhs,lp}/run.R, untouched)
#   scan_fast  : the same statistic bodies with two script-local vectorised
#                primitives swapped in (cumsum-based summed-area table;
#                direct vectorised Gaussian Gram) -- the measured
#                "optimization possibilities" answer. Verified equal to the
#                shipped output before any timing.
#   primitive  : the two O(n^2) build steps timed alone, shipped vs fast,
#                so cost is attributed correctly (the window scan itself is
#                already vectorised; NumMatrix() and gram() are the
#                interpreter-bound passes).
#   pipeline   : end-to-end functional check at max(N_GRID): generate F1
#                curves (168-point grid), project (L2 norm, Ito), scan.
#
# Safeguards: every cell is pre-flighted against an allocation estimate and
# the machine's MemAvailable (skipped, not attempted, if it does not fit),
# every timed call sits in tryCatch, and large objects are freed between
# cells -- a careless N_GRID edit must not freeze the machine.
#
# Results CSV goes to data/ (committed); results/ is gitignored.
# Usage: Rscript reproduction/benchmark_scaling.R   (no null caches needed)

# ---- Inputs (edit these) --------------------------------------------------
QUICK           <- TRUE            # TRUE: sub-minute sanity grid; FALSE: full run
N_GRID_QUICK    <- c(400L, 800L)
N_GRID_FULL     <- c(500L, 1000L, 2000L, 3000L, 4000L, 5000L, 6000L, 7000L)
                                   # upper cells probe where each method hits the
                                   # memory budget; pre-flight skips what won't fit
LP_N_MAX        <- 7000L           # Lp loops over ~n^2/2 windows in R; cells
                                   # above this are skipped and extrapolated O(n^2)
GAMMA           <- 0.25
RKHS_KERNEL     <- "gaussian"
L_GRID          <- 167L            # F1 curves on 168 grid points (weekly-hourly)
SEED            <- 20260822L
TIMING_REPS     <- 3L              # median of; cells slower than the cutoff run once
SLOW_REP_CUTOFF <- 5               # seconds
MEM_BUDGET_MB   <- 4000            # hard per-cell ceiling regardless of free RAM
MEM_SAFETY_FRAC <- 0.5             # also stay under this fraction of MemAvailable
# ---------------------------------------------------------------------------

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- sub(key, "", args[grepl(key, args)])
  if (length(hit)) return(normalizePath(hit[1]))
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) return(normalizePath(ofile))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(script_path)

source(file.path(root, "R", "legacy", "integral_image.R"))
source(file.path(root, "R", "legacy", "weights.R"))
source(file.path(root, "R", "legacy", "kernels.R"))
source(file.path(root, "R", "legacy", "progress.R"))
source(file.path(root, "R", "legacy", "projections.R"))
for (m in c("wilcoxon", "rkhs", "lp")) {
  source(file.path(root, "R", "legacy", paste0("method_", m, ".R")))
}
suppressPackageStartupMessages(library(fdachange)); source(file.path(root, "R", "compat.R"))

N_GRID  <- if (QUICK) N_GRID_QUICK else N_GRID_FULL
# QUICK sanity runs go to gitignored results/; only full runs update the
# committed data/ CSV.
OUT_CSV <- file.path(root, if (QUICK) "results" else "data",
                     "benchmark_limits.csv")

# ---- Fast primitives (script-local; shipped code is not modified) ----------
# Summed-area table by two cumsum passes; same zero-padded (n+1)x(n+1) layout
# as NumMatrix()$sat.
sat_fast <- function(img) {
  cs  <- apply(img, 2L, cumsum)
  cs  <- t(apply(cs, 1L, cumsum))
  sat <- matrix(0, nrow(img) + 1L, ncol(img) + 1L)
  sat[-1L, -1L] <- cs
  sat
}
# Gaussian Gram without the per-pair Vectorize() dispatch of gram().
gram_fast_gaussian <- function(x) exp(-outer(x, x, "-")^2 / 2)

# Shipped statistic bodies re-bound to the fast primitives: measures the full
# scan with the two-line rewrites, with zero duplication of scan logic.
fast_env <- list2env(list(
  NumMatrix    = function(img) list(sat = sat_fast(img)),
  gram         = function(X, kernel) gram_fast_gaussian(X),
  weight_cache = weight_cache), parent = globalenv())
wilcoxon_statistic_fast <- wilcoxon_statistic
environment(wilcoxon_statistic_fast) <- fast_env
rkhs_statistic_fast <- rkhs_statistic
environment(rkhs_statistic_fast) <- fast_env

# ---- Safeguards -------------------------------------------------------------
mem_available_mb <- function() {
  if (!file.exists("/proc/meminfo")) return(NA_real_)
  hit <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(strsplit(trimws(sub("^MemAvailable:", "", hit[1])), "\\s+")[[1]][1]) / 1024
}
effective_budget_mb <- function() {
  avail <- mem_available_mb()
  if (is.finite(avail)) min(MEM_BUDGET_MB, MEM_SAFETY_FRAC * avail) else MEM_BUDGET_MB
}
# Allocation profile per cell kind, in bytes per n^2 (except lp/pipeline).
# Sized to MEASURED peaks at n = 3000, plus margin: Wilcoxon scan ~87 B/n^2,
# fast RKHS ~95, shipped RKHS ~155 (gram()'s per-pair mapply dispatch roughly
# doubles the matrix cost), sat/gram primitives 38-146. The per-kind factors
# decide where each method hits MEM_BUDGET_MB, so keep them honest.
est_cell_mb <- function(kind, n) {
  bytes <- switch(kind,
    scan_wilcoxon = 110 * n^2,
    scan_rkhs     = 160 * n^2,
    scan_fast     = 110 * n^2,
    gram_shipped  = 160 * n^2,
    primitive     = 90 * n^2,
    lp            = 8 * n * 101 * 4 + 4e6,  # ind matrix copies + loop locals
    pipeline      = 8 * (L_GRID + 1) * n * 4 + 4e6,
    stop("unknown kind: ", kind))
  bytes / 1048576
}

# ---- Timing harness ---------------------------------------------------------
RESULTS <- list()
add_row <- function(part, method, n, reps, median_s, peak_mb, value, note) {
  RESULTS[[length(RESULTS) + 1L]] <<- data.frame(
    part = part, method = method, n = n, reps = reps,
    median_s = median_s, peak_mb = peak_mb,
    value = value, note = note, stringsAsFactors = FALSE)
}
time_call <- function(f) {
  gc(reset = TRUE)
  v  <- NULL
  t1 <- system.time(v <- f())[["elapsed"]]
  ts <- t1
  if (t1 <= SLOW_REP_CUTOFF && TIMING_REPS > 1L) {
    for (i in 2:TIMING_REPS) ts <- c(ts, system.time(f())[["elapsed"]])
  }
  g <- gc()
  list(median_s = median(ts), reps = length(ts),
       peak_mb = sum(g[, 6L]), value = v)
}
run_cell <- function(part, method, n, kind, thunk) {
  est <- est_cell_mb(kind, n)
  bud <- effective_budget_mb()
  if (est > bud) {
    cat(sprintf("%-10s %-16s %5d  SKIPPED  needs ~%.0f MB > budget %.0f MB\n",
                part, method, n, est, bud))
    add_row(part, method, n, 0L, NA_real_, NA_real_, NA_real_, "skipped_preflight")
    return(invisible(NULL))
  }
  out <- tryCatch(time_call(thunk), error = identity)
  if (inherits(out, "error")) {
    cat(sprintf("%-10s %-16s %5d  FAILED   %s\n",
                part, method, n, conditionMessage(out)))
    gc()
    add_row(part, method, n, 0L, NA_real_, NA_real_, NA_real_,
            paste0("error: ", substr(conditionMessage(out), 1L, 80L)))
    return(invisible(NULL))
  }
  val  <- suppressWarnings(as.numeric(out$value)[1])
  note <- if (part %in% c("scan", "scan_fast") && !(is.finite(val) && val > 0)) {
    "nonfinite_or_zero_value"
  } else ""
  cat(sprintf("%-10s %-16s %5d  %8s  reps=%d  peak=%7.0f MB  stat=%.4g\n",
              part, method, n, fmt_duration(out$median_s), out$reps,
              out$peak_mb, val))
  add_row(part, method, n, out$reps, out$median_s, out$peak_mb, val, note)
  invisible(out$value)
}

# ---- Correctness gate: fast primitives must reproduce shipped output --------
set.seed(SEED)
x0 <- rnorm(200)
H0 <- outer(x0, x0, function(x, y) sign(x - y))
stopifnot(isTRUE(all.equal(sat_fast(H0), NumMatrix(H0)$sat)))
if (RKHS_KERNEL == "gaussian") {
  stopifnot(isTRUE(all.equal(gram_fast_gaussian(x0), gram(x0, kernel_fn(RKHS_KERNEL)))))
}
stopifnot(isTRUE(all.equal(wilcoxon_statistic_fast(x0, GAMMA),
                           wilcoxon_statistic(x0, GAMMA))))
stopifnot(isTRUE(all.equal(rkhs_statistic_fast(x0, GAMMA, kernel_fn(RKHS_KERNEL)),
                           rkhs_statistic(x0, GAMMA, kernel_fn(RKHS_KERNEL)))))
cat("fast-primitive equality checks passed (n = 200)\n")

cat(sprintf("# %s  QUICK=%s  N_GRID=%s  LP_N_MAX=%d  budget=%.0f MB (avail %.0f MB)\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), QUICK,
            paste(N_GRID, collapse = ","), LP_N_MAX,
            effective_budget_mb(), mem_available_mb()))
cat(sprintf("%-10s %-16s %5s  %8s\n", "part", "method", "n", "median"))
cat(strrep("-", 78), "\n", sep = "")

# ---- Part A + B: scan grid, fast variants, primitive attribution ------------
krn <- kernel_fn(RKHS_KERNEL)
for (n in N_GRID) {
  set.seed(SEED + n)
  X <- rnorm(n)
  run_cell("scan", "wilcoxon", n, "scan_wilcoxon", function() wilcoxon_statistic(X, GAMMA))
  run_cell("scan", "rkhs",     n, "scan_rkhs",     function() rkhs_statistic(X, GAMMA, krn))
  if (n <= LP_N_MAX) {
    run_cell("scan", "lp", n, "lp", function() lp_statistic(X, GAMMA))
  } else {
    cat(sprintf("%-10s %-16s %5d  capped (LP_N_MAX=%d; see extrapolation)\n",
                "scan", "lp", n, LP_N_MAX))
    add_row("scan", "lp", n, 0L, NA_real_, NA_real_, NA_real_, "capped_lp_n_max")
  }
  run_cell("scan_fast", "wilcoxon_fast", n, "scan_fast",
           function() wilcoxon_statistic_fast(X, GAMMA))
  run_cell("scan_fast", "rkhs_fast", n, "scan_fast",
           function() rkhs_statistic_fast(X, GAMMA, krn))
  # The two sat cells share the sign matrix H; do not even build it when the
  # cells cannot run under the budget.
  if (est_cell_mb("primitive", n) <= effective_budget_mb()) {
    H <- outer(X, X, function(x, y) sign(x - y))
    run_cell("primitive", "sat_shipped", n, "primitive", function() NumMatrix(H)$sat)
    run_cell("primitive", "sat_fast",    n, "primitive", function() sat_fast(H))
    rm(H)
  } else {
    for (mm in c("sat_shipped", "sat_fast")) {
      cat(sprintf("%-10s %-16s %5d  SKIPPED  needs ~%.0f MB > budget %.0f MB\n",
                  "primitive", mm, n, est_cell_mb("primitive", n), effective_budget_mb()))
      add_row("primitive", mm, n, 0L, NA_real_, NA_real_, NA_real_, "skipped_preflight")
    }
  }
  run_cell("primitive", "gram_shipped", n, "gram_shipped", function() gram(X, krn))
  run_cell("primitive", "gram_fast",    n, "primitive",    function() gram_fast_gaussian(X))
  rm(X); invisible(gc())
}

# ---- Part C: end-to-end functional pipeline ---------------------------------
# Runs at the largest grid n whose Wilcoxon scan fits the memory budget, so the
# end-to-end check is a measurement, not a pre-flight skip.
n_fit  <- N_GRID[vapply(N_GRID, function(n)
  est_cell_mb("scan_wilcoxon", n) <= effective_budget_mb(), logical(1))]
n_pipe <- if (length(n_fit)) max(n_fit) else min(N_GRID)
cat(strrep("-", 78), "\n", sep = "")
cat(sprintf("pipeline: n = %d functional observations on %d grid points%s\n",
            n_pipe, L_GRID + 1L,
            if (n_pipe < max(N_GRID))
              sprintf(" (largest grid n fitting the %.0f MB budget)",
                      effective_budget_mb()) else ""))
set.seed(SEED)
Xf <- NULL
run_cell("pipeline", "generate_F1", n_pipe, "pipeline",
         function() { Xf <<- generate_F1(n_pipe, a = 0, l = L_GRID); ncol(Xf) })
if (!is.null(Xf)) {
  x_l2 <- NULL
  run_cell("pipeline", "project_l2",  n_pipe, "pipeline",
           function() { x_l2 <<- project_l2(Xf); length(x_l2) })
  run_cell("pipeline", "project_ito", n_pipe, "pipeline",
           function() length(project_ito(Xf)))
  if (!is.null(x_l2)) {
    run_cell("pipeline", "wilcoxon_on_l2",  n_pipe, "scan_wilcoxon",
             function() wilcoxon_statistic(x_l2, GAMMA))
    run_cell("pipeline", "rkhs_fast_on_l2", n_pipe, "scan_fast",
             function() rkhs_statistic_fast(x_l2, GAMMA, krn))
  }
  rm(Xf); invisible(gc())
}

# ---- Summary: O(n^2) sanity + Lp extrapolation -------------------------------
df <- do.call(rbind, RESULTS)
cat(strrep("-", 78), "\n", sep = "")
w <- df[df$part == "scan" & df$method == "wilcoxon" & is.finite(df$median_s), ]
if (nrow(w) >= 2L) {
  for (i in 2:nrow(w)) {
    cat(sprintf("O(n^2) sanity, wilcoxon %d->%d: time x%.1f vs n^2 x%.1f\n",
                w$n[i - 1], w$n[i], w$median_s[i] / w$median_s[i - 1],
                (w$n[i] / w$n[i - 1])^2))
  }
}
lp <- df[df$part == "scan" & df$method == "lp" & is.finite(df$median_s), ]
if (nrow(lp) && max(N_GRID) > max(lp$n)) {
  t_cap <- lp$median_s[which.max(lp$n)]
  n_cap <- max(lp$n)
  cat(sprintf("lp O(n^2) extrapolation: %s at n=%d -> ~%s at n=%d\n",
              fmt_duration(t_cap), n_cap,
              fmt_duration(t_cap * (max(N_GRID) / n_cap)^2), max(N_GRID)))
}
for (m in c("wilcoxon", "rkhs", "lp")) {
  ran <- df[df$part == "scan" & df$method == m & is.finite(df$median_s), ]
  if (nrow(ran)) {
    i <- which.max(ran$n)
    cat(sprintf("largest n under the %.0f MB budget, %-8s: %d (%s, peak %.0f MB)\n",
                MEM_BUDGET_MB, m, ran$n[i], fmt_duration(ran$median_s[i]),
                ran$peak_mb[i]))
  }
}

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df, OUT_CSV, row.names = FALSE)
cat("wrote ", OUT_CSV, "\n", sep = "")
