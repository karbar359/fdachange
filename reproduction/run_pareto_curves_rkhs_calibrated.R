#!/usr/bin/env Rscript
# Re-run RKHS for F1/F2/F3 with proper baseline alignment.
#
# The asymptotic n=2000 cache (.raw/SD_h_n2000_dist_t.csv) is calibrated
# against specific baseline distributions:
#   gaussian (no suffix)     <- N(0, 1)
#   min .runif               <- U(0, 1)
#   min .rnorm.sq            <- N(0, 1)^2  (chi-squared 1 d.f.)
#   2_max .runif             <- U(0, 1)
#
# To apply these caches to projected functional samples we mirror the
# T1/T2 pipeline of the litgrid application (sec_framework.tex):
#   T1: F_n (rank/n) -> U(0, 1)
#   T2g: Q_norm(F_n(.)) -> N(0, 1)        (for the gaussian kernel)
#
# For each (projection, kernel) pair this script picks the matching
# baseline column so the rejection decision is calibrated:
#   ito        -> Gaussian / N(0, 1)
#   l2         -> min / U(0, 1)
#   c_max_var  -> Gaussian / N(0, 1)
#
# Outputs land in `results/rkhs_<scenario>_<tag>_<projection>.csv`,
# overwriting the previous (under-calibrated) files for these cells.

GAMMA_GRID <- c(0, 0.25, 0.45)
# Per-scenario alternative grids (matching the dissertation prose):
#   F1: a in {1.5, 2, 2.5, 3, 5, 10}  (a = 1.5 is the size check, shape unchanged)
#   F2: a in {0.5, 1, 1.5, 2, 5}      (a ~ 1.5 is the implicit size check;
#                                       centred Pareto(1, 1.5) has variance ~ N(0, 1))
#   F3: a in {0.5, 1, 1.5, 2, 5}      (same convention as F2)
A_GRID_BY_SCENARIO <- list(
  F1 = c(1.5, 2, 2.5, 3, 5, 10),
  F2 = c(0.5, 1, 1.5, 2, 5),
  F3 = c(0.5, 1, 1.5, 2, 2.5, 5)
)
N          <- 200L
N_REP      <- 100L
ALPHA      <- 0.05
NCORES     <- 12L
SEED       <- 42L

# Centred-window epidemic geometry, identical to driver.R defaults at n=200.
EC_GEOMETRY <- (function() {
  n <- N
  l_ints <- c(floor(log(n)) + 1L, n %/% 5L)
  lapply(l_ints, function(l) {
    k_star <- (n %/% 2L) - (l %/% 2L)
    c(k_star / n, (k_star + l) / n)
  })
})()

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
res_dir <- file.path(root, "results")

source(file.path(root, "R", "legacy", "integral_image.R"))
source(file.path(root, "R", "legacy", "weights.R"))
source(file.path(root, "R", "legacy", "kernels.R"))
source(file.path(root, "R", "legacy", "null_cache.R"))
source(file.path(root, "R", "legacy", "projections.R"))
source(file.path(root, "R", "legacy", "method_rkhs.R"))

suppressPackageStartupMessages(library(parallel))

# (projection, kernel, baseline_tag, transform_fn) per cell.
# `transform_fn(x)` maps the projected scalar series to whatever marginal
# the cached null was calibrated on.
T1_uniform <- function(x) rank(x, ties.method = "average") / length(x)
T2_normal  <- function(x) qnorm(pmin(pmax(T1_uniform(x), 1e-9), 1 - 1e-9))

CELLS <- list(
  list(projection = "ito",       kernel = "gaussian", baseline = "",      transform = T2_normal),
  list(projection = "l2",        kernel = "min",      baseline = "runif", transform = T1_uniform),
  list(projection = "c_max_var", kernel = "gaussian", baseline = "",      transform = T2_normal)
)

run_cell <- function(scenario, cell, a_grid) {
  # scenario generators come from R/compat.R (fdachange package)
  generate <- get(paste0("generate_", scenario))

  proj_fn <- projection_fn(cell$projection)
  cat(sprintf("[rkhs/%s proj=%s kernel=%s baseline=%s a_grid=(%s)] starting\n",
              scenario, cell$projection, cell$kernel,
              if (nzchar(cell$baseline)) cell$baseline else "(N(0,1))",
              paste(a_grid, collapse = ", ")))

  rows <- list()
  for (gi in seq_along(GAMMA_GRID)) {
    gamma <- GAMMA_GRID[gi]
    null_draws <- rkhs_null_distribution(
      N, gamma, kernel_name = cell$kernel, baseline = cell$baseline)
    cv <- critical_value(null_draws, ALPHA, m.d = 1L)
    for (eci in seq_along(EC_GEOMETRY)) {
      ec <- EC_GEOMETRY[[eci]]
      ec_length <- ec[2] - ec[1]
      for (a in a_grid) {
        set.seed(SEED + 1000L * gi + 100L * eci + as.integer(round(a * 4)))
        powers <- mclapply(seq_len(N_REP), function(rep) {
          X <- generate(N, a = a, ec.interval = ec)
          x_proj <- proj_fn(X)
          x_in <- cell$transform(x_proj)
          stat <- rkhs_statistic(x_in, gamma, kernel = kernel_fn(cell$kernel))
          stat > cv
        }, mc.cores = NCORES)
        power <- mean(unlist(powers))
        rows[[length(rows) + 1L]] <- data.frame(
          method = "rkhs", scenario = scenario,
          tag = paste0(cell$kernel,
                       if (nzchar(cell$baseline)) paste0("_", cell$baseline) else "",
                       "_", cell$projection),
          n = N, m = 1L, projection = cell$projection,
          gamma = gamma, alpha = ALPHA, N_rep = N_REP,
          ec_start = ec[1], ec_end = ec[2], ec_length = ec_length,
          a = a, power = power, critical_value = cv)
      }
    }
  }
  out <- do.call(rbind, rows)
  out_csv <- file.path(res_dir, sprintf("rkhs_%s_%s.csv", scenario, out$tag[1]))
  write.csv(out, out_csv, row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows)\n", out_csv, nrow(out)))
}

# Clean previous (under-calibrated) RKHS files for F1/F2/F3 so the
# aggregator only sees the new ones.
old <- list.files(res_dir, pattern = "^rkhs_F[123]_.*\\.csv$", full.names = TRUE)
if (length(old)) {
  cat(sprintf("Removing %d stale rkhs F1/F2/F3 csv(s)\n", length(old)))
  file.remove(old)
}

for (scenario in names(A_GRID_BY_SCENARIO)) {
  a_grid <- A_GRID_BY_SCENARIO[[scenario]]
  for (cell in CELLS) {
    run_cell(scenario, cell, a_grid)
  }
}

cat(sprintf("\n=== %s : 9 RKHS cells (3 scenarios x 3 proj/kernel) complete ===\n",
            format(Sys.time(), "%H:%M:%S")))
