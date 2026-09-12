#!/usr/bin/env Rscript
# How much does the changed-segment conclusion depend on the curve
# reconstruction (smoothing) step?
#
# The dissertation's English summary claims the conclusions are "qualitatively
# stable under reasonable changes of preprocessing choices (e.g., smoothing
# level, basis truncation, and projection dimension)"
# (dissertation/summary/simulation_results.tex:34), but no such experiment
# exists anywhere in the repository. This script is that experiment.
#
# Two arms:
#   application -- the Litgrid weekly balancing-price curves, sweeping the
#                  reconstruction across the three reported projections and the
#                  gamma grid, reporting statistic, critical value, rejection
#                  and detected window.
#   simulation  -- scenario F1 at the null and at two alternatives, reporting
#                  empirical size and power under the same sweep.
#
# Reconstruction methods (shared/reconstruction.R): `none` is the raw sampling
# grid, `fs_interp` is Faber-Schauder truncation (what the thesis text
# describes), `fs_ls` is least squares onto the same polygonal space (what the
# original code does), `bspline_pen` is a genuinely smooth cubic contrast.
#
# Uses `wilcoxon_scan_corrected` (2026-09 audit: the printed statistic with
# the weight applied in the right direction; the earlier `_exact` variant
# fixed only the summed-area off-by-one and still carried the inverted
# weight, so pre-audit copies of the committed CSV understate gamma > 0).
#
# Usage:  Rscript reproduction/smoothing_sensitivity.R

# ---- Inputs (edit these) --------------------------------------------------
QUICK        <- TRUE          # TRUE keeps the whole run under a minute
DEDUP        <- TRUE          # deduplicate the overlapping raw exports
PROJECTIONS  <- c("l2", "ito", "c_max_var")
GAMMAS       <- c(0.05, 0.25, 0.45)
K_GRID       <- c(3L, 4L, 5L, 6L, 7L)
N_DRAWS      <- if (QUICK) 3L else 10L   # redraws for random projections
ALPHA        <- 0.05
KERNEL       <- "h_0"
SEED         <- 3L
# simulation arm
SIM_N        <- if (QUICK) 100L else 200L
SIM_NREP     <- if (QUICK) 20L  else 200L
SIM_A        <- c(0, 2, 5)     # 0 = size, others = power
SIM_EC       <- c(0.4, 0.7)
# Litgrid trim, matching the dissertation
DROP_WEEK    <- 202149L   # outlier week, BY LABEL (fixed pipeline moves indices)
SPAN_END     <- 202433L   # last week the printed trim covered
OUT_CSV      <- "data/smoothing_sensitivity.csv"
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(fda); library(fda.usc) })

root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) dirname(dirname(normalizePath(hit[1]))) else getwd()
})
source(file.path(root, "R", "legacy", "integral_image.R"))
source(file.path(root, "R", "legacy", "weights.R"))
source(file.path(root, "R", "legacy", "null_cache.R"))
source(file.path(root, "R", "legacy", "projections.R"))
source(file.path(root, "R", "legacy", "reconstruction.R"))
source(file.path(root, "R", "legacy", "progress.R"))
source(file.path(root, "R", "legacy", "method_wilcoxon.R"))
source(file.path(root, "applications", "litgrid", "data.R"))
suppressPackageStartupMessages(library(fdachange)); source(file.path(root, "R", "compat.R"))

# The reconstruction sweep. `vertices` evaluates on the 2^K + 1 polygon
# vertices (thesis-faithful); `input` keeps the original sampling grid, which
# isolates the effect of the evaluation grid from that of the truncation level.
RECON_SPECS <- c(
  list(list(method = "none", K = NA_integer_, grid = "input")),
  unlist(lapply(K_GRID, function(k) list(
    list(method = "fs_interp",   K = k, grid = "vertices"),
    list(method = "fs_ls",       K = k, grid = "vertices"),
    list(method = "fs_interp",   K = k, grid = "input"),
    list(method = "fs_ls",       K = k, grid = "input"))), recursive = FALSE),
  list(list(method = "bspline_pen", K = 4L, grid = "input"),
       list(method = "bspline_pen", K = 6L, grid = "input")))

spec_label <- function(s) recon_label(s$method, s$K, s$grid)

rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

is_random <- function(projection) projection %in% c("ito", "c_max_var", "wiener")

# ============================ application arm ==============================
litgrid_dir <- Sys.getenv("LITGRID_RAW_DIR", unset = "")
if (!nzchar(litgrid_dir)) {
  litgrid_dir <- normalizePath(file.path(root, "..", "data-raw", "litgrid-raw"),
                               mustWork = FALSE)
}
if (!dir.exists(litgrid_dir)) {
  stop("Litgrid raw CSVs not found at ", litgrid_dir,
       "; run Rscript applications/litgrid/download_litgrid.R or set LITGRID_RAW_DIR.")
}
litgrid_dir <- litgrid_dir

cat("== application arm: Litgrid weekly balancing prices ==\n")
f_matrix <- suppressWarnings(
  load_litgrid_balancing_prices(litgrid_dir, dedup_exports = DEDUP))
lab <- as.integer(rownames(f_matrix))
f_matrix <- f_matrix[lab <= SPAN_END & lab != DROP_WEEK, , drop = FALSE]
X_app <- t(f_matrix)
week_labels <- rownames(f_matrix)
n_app <- ncol(X_app)
cat(sprintf("   %d weeks x %d hours, dedup = %s\n", n_app, nrow(X_app), DEDUP))

cv_app <- vapply(GAMMAS, function(g)
  critical_value(wilcoxon_null_distribution(n_app, g), ALPHA, 1L), numeric(1))
names(cv_app) <- as.character(GAMMAS)

total_app <- length(RECON_SPECS) * length(PROJECTIONS)
done <- 0L
for (spec in RECON_SPECS) {
  Xr <- reconstruct(X_app, spec$method, if (is.na(spec$K)) 4L else spec$K,
                    grid = spec$grid)$X
  for (proj in PROJECTIONS) {
    done <- done + 1L
    cat(sprintf("   [%2d/%2d] %-22s x %-10s\n", done, total_app, spec_label(spec), proj))
    utils::flush.console()
    project <- projection_fn(proj)
    n_draw  <- if (is_random(proj)) N_DRAWS else 1L
    for (g in GAMMAS) {
      set.seed(SEED)
      st <- numeric(n_draw); lo <- numeric(n_draw); hi <- numeric(n_draw)
      for (b in seq_len(n_draw)) {
        v  <- project(Xr)
        sc <- wilcoxon_scan_corrected(v, g, kernel = KERNEL)
        st[b] <- sc$stat; lo[b] <- sc$interval[1]; hi[b] <- sc$interval[2]
      }
      add_row(arm = "application", recon = spec_label(spec),
              method = spec$method, K = spec$K, grid = spec$grid,
              projection = proj, gamma = g, n = n_app,
              stat = mean(st), stat_sd = if (n_draw > 1) sd(st) else NA_real_,
              cv = cv_app[[as.character(g)]],
              reject_rate = mean(st > cv_app[[as.character(g)]]),
              lo = mean(lo), hi = mean(hi),
              lo_sd = if (n_draw > 1) sd(lo) else NA_real_,
              hi_sd = if (n_draw > 1) sd(hi) else NA_real_,
              week_lo = week_labels[max(1L, round(mean(lo) * n_app))],
              week_hi = week_labels[min(n_app, round(mean(hi) * n_app))],
              a = NA_real_, nrep = n_draw)
    }
  }
}

# ============================ simulation arm ===============================
cat(sprintf("== simulation arm: F1, n = %d, nrep = %d ==\n", SIM_N, SIM_NREP))
cv_sim <- vapply(GAMMAS, function(g)
  critical_value(wilcoxon_null_distribution(SIM_N, g), ALPHA, 1L), numeric(1))
names(cv_sim) <- as.character(GAMMAS)

sim_specs <- Filter(function(s) s$grid == "vertices" || s$method %in% c("none", "bspline_pen"),
                    RECON_SPECS)
total_sim <- length(sim_specs) * length(SIM_A)
done <- 0L
for (spec in sim_specs) {
  for (a in SIM_A) {
    done <- done + 1L
    cat(sprintf("   [%2d/%2d] %-22s a = %-4g\n", done, total_sim, spec_label(spec), a))
    utils::flush.console()
    set.seed(SEED)
    hits <- matrix(0, nrow = SIM_NREP, ncol = length(GAMMAS))
    for (r in seq_len(SIM_NREP)) {
      Xs <- generate_F1(SIM_N, a = a, ec.interval = SIM_EC)
      Xr <- reconstruct(Xs, spec$method, if (is.na(spec$K)) 4L else spec$K,
                        grid = spec$grid)$X
      v  <- project_l2(Xr)
      for (gi in seq_along(GAMMAS)) {
        hits[r, gi] <- wilcoxon_scan_corrected(v, GAMMAS[gi], kernel = KERNEL)$stat >
                       cv_sim[[as.character(GAMMAS[gi])]]
      }
    }
    for (gi in seq_along(GAMMAS)) {
      add_row(arm = if (a == 0) "simulation_size" else "simulation_power",
              recon = spec_label(spec), method = spec$method, K = spec$K,
              grid = spec$grid, projection = "l2", gamma = GAMMAS[gi],
              n = SIM_N, stat = NA_real_, stat_sd = NA_real_,
              cv = cv_sim[[as.character(GAMMAS[gi])]],
              reject_rate = mean(hits[, gi]),
              lo = NA_real_, hi = NA_real_, lo_sd = NA_real_, hi_sd = NA_real_,
              week_lo = NA_character_, week_hi = NA_character_,
              a = a, nrep = SIM_NREP)
    }
  }
}

out <- do.call(rbind, rows)
out_path <- if (QUICK) {
  file.path(root, "results", basename(OUT_CSV))
} else {
  file.path(root, OUT_CSV)
}
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_path, row.names = FALSE)
cat(sprintf("\nwrote %d rows to %s%s\n", nrow(out), out_path,
            if (QUICK) "  (QUICK run; set QUICK <- FALSE for the committed table)" else ""))
