#!/usr/bin/env Rscript
# Stability of random projections in the Wilcoxon-type test, ported onto the
# fdachange package (scenarios and the Ito projection) with the statistic
# variant selected by WILCOXON_SCAN: legacy/exact from R/legacy/, corrected
# from the package.
#
# Following the simulations chapter, we *fix* a sample X = (X_1, ..., X_n)
# of n functional observations and ask how the test statistic T_n(omega)
# varies under independent draws of the Ito-integral projector W. For each
# fixed sample we report the empirical rejection rate
#
#     R_hat_N = (1/N) sum_i 1{T_n(omega_i) > c_alpha}.
#
# This is structurally different from `driver.R`, which redraws both X and W
# per replicate and reports the joint power; here we redraw only W on top
# of a frozen X (and optionally repeat for K independently drawn fixed X
# samples to also see between-sample variability).
#
# Output:
#   - `results/stability_raw.csv`     -- one row per (cell, draw, gamma)
#   - `results/stability_summary.csv` -- one row per (cell, gamma)

# ---- Inputs (edit these, or pre-set in .GlobalEnv from a wrapper) ---------
.def <- function(name, value) {
  if (!exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    assign(name, value, envir = .GlobalEnv)
  }
}

.def("N",          400L)
.def("K",          100L)                      # independent fixed X samples per (scenario, a)
.def("N_PROJ",     100L)                      # independent Ito draws per fixed X

.def("SCENARIOS",  c("F1", "F2", "F3", "F4"))

# Scenario-specific a-grid: at least one size check + 1-2 alternatives.
.def("A_BY_SCENARIO", list(
  F1 = c(2.5, 5, 10),
  F2 = c(1.5, 2, 5),
  F3 = c(2.5, 5),
  F4 = c(0, 1, 4)))

.def("GAMMA_GRID", c(0, 0.25, 0.45))
.def("ALPHA",      0.05)
.def("M",          1L)
.def("EC_FRACTION", c(0.4, 0.6))              # epidemic window I^* = (0.4, 0.6) * n

.def("BETA",            0)
# Norm-based Wilcoxon kernel h_0(x,y)=sign(|x|-|y|); the functional scenarios
# shift the scale/tail of the in-window marginal, which h_w is invariant to.
.def("WILCOXON_KERNEL", "h_0")
# "legacy" is the scan that produced the printed stability table; "exact"
# fixes only the summed-area off-by-one; "corrected" is the printed statistic
# (fdachange package). Set OUT_RAW / OUT_SUMMARY when switching so the
# families do not share an output file.
.def("WILCOXON_SCAN", "legacy")               # legacy | exact | corrected

.def("SEED",       42L)
.def("NCORES",     8L)
.def("PROGRESS_EVERY", 20L)

.def("OUT_RAW",     NULL)                     # NULL -> results/stability_raw.csv
.def("OUT_SUMMARY", NULL)                     # NULL -> results/stability_summary.csv
.def("HEARTBEAT",   NULL)                     # NULL -> results/stability_heartbeat.log

# Smoke override: STABILITY_SMOKE=1 in env -> shrink everything to <30s.
if (Sys.getenv("STABILITY_SMOKE") == "1") {
  N <- 80L; K <- 1L; N_PROJ <- 20L
  SCENARIOS <- "F1"; A_BY_SCENARIO <- list(F1 = c(2.5, 10))
  GAMMA_GRID <- c(0); NCORES <- 1L; PROGRESS_EVERY <- 4L
  OUT_RAW     <- file.path("results", "stability_raw_smoke.csv")
  OUT_SUMMARY <- file.path("results", "stability_summary_smoke.csv")
}

rm(.def, envir = .GlobalEnv)
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

suppressPackageStartupMessages(library(fdachange))
source(file.path(root, "R", "compat.R"))
source(file.path(root, "R", "legacy_env.R"))
legacy <- legacy_env(root)

if (!WILCOXON_SCAN %in% c("legacy", "exact", "corrected")) {
  stop("WILCOXON_SCAN must be 'legacy', 'exact' or 'corrected', got '",
       WILCOXON_SCAN, "'")
}
# The corrected batched-gamma call is the package's vectorised-gamma
# statistic; the pre-audit families come from the vendored code.
.wilcoxon_gammas_fn <- switch(WILCOXON_SCAN,
  legacy    = legacy$wilcoxon_statistic_gammas,
  exact     = legacy$wilcoxon_statistic_gammas_exact,
  corrected = fdachange::wilcoxon_statistic)

if (is.null(OUT_RAW))     OUT_RAW     <- file.path(root, "results", "stability_raw.csv")
if (is.null(OUT_SUMMARY)) OUT_SUMMARY <- file.path(root, "results", "stability_summary.csv")
if (is.null(HEARTBEAT))   HEARTBEAT   <- file.path(root, "results", "stability_heartbeat.log")
dir.create(dirname(OUT_RAW),     showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_SUMMARY), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(HEARTBEAT),   showWarnings = FALSE, recursive = TRUE)
writeLines(sprintf("[%s] start  n=%d  K=%d  N_PROJ=%d  scenarios=%s",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   N, K, N_PROJ, paste(SCENARIOS, collapse = ",")),
           HEARTBEAT)

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

# Pre-compute critical values for each gamma. Asymptotic Wilcoxon null is
# n-free, so the cv depends only on (alpha, gamma, m.d=M). FAITHFUL PORT
# NOTE: like the original, the cv always comes from the cached null with
# exact_gamma0 = FALSE, i.e. even corrected-scan stability runs used the
# legacy gamma = 0 cache cell -- preserved so the corrected summary CSV
# reproduces bit for bit.
cv_by_gamma <- vapply(GAMMA_GRID,
                      function(g) legacy$wilcoxon_critical_value(
                        N, g, beta = BETA, alpha = ALPHA, m.d = M),
                      numeric(1))
names(cv_by_gamma) <- as.character(GAMMA_GRID)
cat(sprintf("[critical values @ alpha=%g, m=%d, beta=%g, kernel=%s]\n",
            ALPHA, M, BETA, WILCOXON_KERNEL))
for (g in GAMMA_GRID) {
  cat(sprintf("  gamma=%.2f  cv=%.4f\n", g, cv_by_gamma[as.character(g)]))
}

# Cell list: (scenario, a, k_sample). Each cell fixes one X and runs N_PROJ
# Ito draws on it.
cells <- list()
for (sc in SCENARIOS) {
  a_grid <- A_BY_SCENARIO[[sc]]
  if (is.null(a_grid)) stop("No A_BY_SCENARIO entry for ", sc)
  for (a in a_grid) {
    for (k in seq_len(K)) {
      cells[[length(cells) + 1L]] <- list(scenario = sc, a = a, k = k)
    }
  }
}
n_cells <- length(cells)
cat(sprintf("[cells] %d cells (scenarios=%s, K=%d, n=%d, n_proj=%d)\n",
            n_cells, paste(SCENARIOS, collapse = ","), K, N, N_PROJ))

run_cell <- function(cell_idx) {
  cell <- cells[[cell_idx]]
  sc <- cell$scenario; a <- cell$a; k <- cell$k
  gen <- get(paste0("generate_", sc))
  set.seed(SEED + cell_idx)                  # reproducible per cell
  X <- gen(N, a, EC_FRACTION)                # FIXED across the inner loop
  T_n <- matrix(NA_real_, nrow = N_PROJ, ncol = length(GAMMA_GRID))
  colnames(T_n) <- as.character(GAMMA_GRID)
  t0 <- Sys.time()
  for (i in seq_len(N_PROJ)) {
    Xp <- project_ito(X)                     # fresh W
    # Batched call: the gamma-independent work (H, prefix sums, sigma) is
    # shared across all gammas in one call.
    T_n[i, ] <- .wilcoxon_gammas_fn(
      Xp, GAMMA_GRID, beta = BETA, kernel = WILCOXON_KERNEL)
    if (NCORES == 1L && (i %% PROGRESS_EVERY == 0L || i == N_PROJ)) {
      el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta <- el * (N_PROJ - i) / max(i, 1L)
      cat(sprintf("[%s a=%g k=%d] i=%d/%d  T_n[g=%.2f]=%.3f  el=%s  eta=%s\n",
                  sc, a, k, i, N_PROJ, GAMMA_GRID[1L], T_n[i, 1L],
                  fmt_duration(el), fmt_duration(eta)))
    }
  }
  cell_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # Append a heartbeat line (O_APPEND writes are atomic for short lines).
  cat(sprintf("[%s] %d/%d  %s a=%g k=%d  done in %.1fs\n",
              format(Sys.time(), "%H:%M:%S"),
              cell_idx, n_cells, sc, a, k, cell_s),
      file = HEARTBEAT, append = TRUE)
  list(cell_idx = cell_idx, scenario = sc, a = a, k = k,
       T_n = T_n, cell_s = cell_s)
}

t0_total <- Sys.time()
results <- if (NCORES > 1L) {
  # mc.preschedule = TRUE: divide cells into NCORES batches, fork once per
  # batch. With K=100 the ~1100 cells are homogeneous, so batching wins.
  parallel::mclapply(seq_len(n_cells), run_cell,
                     mc.cores = NCORES, mc.set.seed = TRUE,
                     mc.preschedule = TRUE)
} else {
  lapply(seq_len(n_cells), run_cell)
}
if (any(vapply(results, inherits, logical(1), "try-error"))) {
  stop("a stability cell errored; rerun with NCORES=1 for a traceback")
}

# Assemble raw long-format CSV (one row per (cell, draw, gamma)) and a
# summary CSV (one row per (cell, gamma)).
raw_rows <- vector("list", n_cells)
sum_rows <- vector("list", n_cells)
for (res in results) {
  T_n   <- res$T_n
  cv_v  <- cv_by_gamma[colnames(T_n)]
  raw_chunks <- lapply(seq_along(GAMMA_GRID), function(gi) {
    g <- GAMMA_GRID[gi]
    data.frame(
      scenario = res$scenario, a = res$a, k_sample = res$k,
      i_draw = seq_len(N_PROJ), gamma = g, T_n = T_n[, gi],
      stringsAsFactors = FALSE)
  })
  raw_rows[[res$cell_idx]] <- do.call(rbind, raw_chunks)

  sum_chunks <- lapply(seq_along(GAMMA_GRID), function(gi) {
    g  <- GAMMA_GRID[gi]
    Tg <- T_n[, gi]
    cv <- cv_v[gi]
    R_hat <- mean(Tg > cv)
    data.frame(
      scenario = res$scenario, a = res$a, k_sample = res$k, gamma = g,
      mean   = mean(Tg), sd = sd(Tg),
      median = median(Tg),
      q25    = quantile(Tg, 0.25, names = FALSE),
      q75    = quantile(Tg, 0.75, names = FALSE),
      R_hat  = R_hat,
      R_hat_se = sqrt(R_hat * (1 - R_hat) / N_PROJ),
      cv     = cv, n = N, n_proj = N_PROJ,
      stringsAsFactors = FALSE)
  })
  sum_rows[[res$cell_idx]] <- do.call(rbind, sum_chunks)

  R_at_g0 <- mean(T_n[, 1L] > cv_v[1L])
  cat(sprintf("[%s a=%g k=%d] cell done in %s  (R_hat[g=%.2f]=%.3f)\n",
              res$scenario, res$a, res$k, fmt_duration(res$cell_s),
              GAMMA_GRID[1L], R_at_g0))
}
total_s <- as.numeric(difftime(Sys.time(), t0_total, units = "secs"))
cat(sprintf("[total] %s\n", fmt_duration(total_s)))

raw_df <- do.call(rbind, raw_rows)
sum_df <- do.call(rbind, sum_rows)

write.csv(raw_df, OUT_RAW, row.names = FALSE)
write.csv(sum_df, OUT_SUMMARY, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", OUT_RAW,     nrow(raw_df)))
cat(sprintf("wrote %s (%d rows)\n", OUT_SUMMARY, nrow(sum_df)))
