#!/usr/bin/env Rscript
# RKHS F1/F2 calibrated rerun at n = 400, N_rep = 100.
#
# Copy-adapt of run_pareto_curves_rkhs_calibrated_n400.R (the pipeline that
# produced the printed thesis tables _pareto_F{1,2}_rkhs_*.tex at N_rep = 20):
#   * scenarios F1/F2 only (F3 is already N_rep = 100 in the thesis tables)
#   * N_REP = 100
#   * outputs carry a thesis_ prefix and _n400_nrep100 tag so
#     aggregate_pareto_curves_tex.R (glob ^(lp|wilcoxon|rkhs)_F[1-3]) can
#     never pick them up when regenerating thesis tables
#   * idempotent: cells whose output CSV already exists are skipped
#   * SMOKE <- TRUE gives a sub-minute end-to-end check (N_REP = 5,
#     first/last a only, smoke_-prefixed outputs)
#
# The driver.R RKHS path is a DIFFERENT statistic pipeline (no rank
# calibration); do not use it for these cells.

# ---- Constants (edit here) ---------------------------------------------
SMOKE      <- identical(Sys.getenv("SMOKE"), "1")  # SMOKE=1 Rscript ... for a sub-minute check
GAMMA_GRID <- c(0, 0.25, 0.45)
A_GRID_BY_SCENARIO <- list(
  F1 = c(1.5, 2, 2.5, 3, 5, 10),
  F2 = c(0.5, 1, 1.5, 2, 5)
)
N          <- 400L
N_REP      <- 100L
ALPHA      <- 0.05
NCORES     <- 16L
SEED       <- 42L
OUT_PREFIX <- "thesis_rkhs_"

if (SMOKE) {
  N_REP <- 5L
  A_GRID_BY_SCENARIO <- lapply(A_GRID_BY_SCENARIO, function(a) a[c(1L, length(a))])
  OUT_PREFIX <- "smoke_thesis_rkhs_"
}

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

T1_uniform <- function(x) rank(x, ties.method = "average") / length(x)
T2_normal  <- function(x) qnorm(pmin(pmax(T1_uniform(x), 1e-9), 1 - 1e-9))

CELLS <- list(
  list(projection = "ito",       kernel = "gaussian", baseline = "",      transform = T2_normal),
  list(projection = "l2",        kernel = "min",      baseline = "runif", transform = T1_uniform),
  list(projection = "c_max_var", kernel = "gaussian", baseline = "",      transform = T2_normal)
)

run_cell <- function(scenario, cell, a_grid) {
  tag <- paste0(cell$kernel,
                if (nzchar(cell$baseline)) paste0("_", cell$baseline) else "",
                "_", cell$projection)
  out_csv <- file.path(res_dir, sprintf("%s%s_n400_nrep100_%s.csv", OUT_PREFIX, scenario, tag))
  if (file.exists(out_csv) && length(readLines(out_csv, n = 2)) > 1) {
    cat(sprintf("SKIP %s (exists)\n", basename(out_csv)))
    return(invisible(NULL))
  }
  # scenario generators come from R/compat.R (fdachange package)
  generate <- get(paste0("generate_", scenario))

  proj_fn <- projection_fn(cell$projection)
  n_cells <- length(GAMMA_GRID) * length(EC_GEOMETRY) * length(a_grid)
  cat(sprintf("[rkhs/%s n=%d N_rep=%d proj=%s kernel=%s baseline=%s a_grid=(%s)] starting (%d cells)\n",
              scenario, N, N_REP, cell$projection, cell$kernel,
              if (nzchar(cell$baseline)) cell$baseline else "(N(0,1))",
              paste(a_grid, collapse = ", "), n_cells))

  done <- 0L
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
          method = "rkhs", scenario = scenario, tag = tag,
          n = N, m = 1L, projection = cell$projection,
          gamma = gamma, alpha = ALPHA, N_rep = N_REP,
          ec_start = ec[1], ec_end = ec[2], ec_length = ec_length,
          a = a, power = power, critical_value = cv)
        done <- done + 1L
        cat(sprintf("  [%s %s] %d/%d gamma=%.2f ec=%.4f a=%g power=%.2f (%s)\n",
                    scenario, tag, done, n_cells, gamma, ec_length, a, power,
                    format(Sys.time(), "%H:%M:%S")))
      }
    }
  }
  out <- do.call(rbind, rows)
  write.csv(out, out_csv, row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows)\n", out_csv, nrow(out)))
}

for (scenario in names(A_GRID_BY_SCENARIO)) {
  a_grid <- A_GRID_BY_SCENARIO[[scenario]]
  for (cell in CELLS) {
    run_cell(scenario, cell, a_grid)
  }
}

cat(sprintf("\n=== %s : thesis RKHS F1/F2 n=400 N_rep=%d sweep complete ===\n",
            format(Sys.time(), "%H:%M:%S"), N_REP))
