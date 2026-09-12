#!/usr/bin/env bash
# Validation re-run for the F1<->F2 identity swap. Runs F1 and F2 only
# (F3 DGP unchanged so its results are reusable). N_REP=20 keeps each
# driver invocation sub-minute; absolute power numbers are MC-noisy at
# this rep count but the qualitative signal (size near alpha at the null,
# rising power at large a) is enough to validate the swap end-to-end.
#
#   F1 (Normal-vs-Pareto, cosine):  a in {1.5, 2, 2.5, 3, 5, 10}
#   F2 (Pareto-vs-Pareto, cosine):  a in {0.5, 1, 1.5, 2, 5}
#
# Wilcoxon: 3 projections x m in {1,2,4} x 2 scenarios = 18 cells
# L_p     : 3 projections x p=2     x 2 scenarios = 6 cells
# RKHS    : 3 projection-kernel pairs x 2 scenarios = 6 cells
# Total: 30 cells. Expect ~15-30 min on 12 cores.

set -euo pipefail

# ---- Inputs (edit these to retune) -----------------------------------------
N=400
N_REP=20
NCORES=12
A_F1='c(1.5, 2, 2.5, 3, 5, 10)'
A_F2='c(0.5, 1, 1.5, 2, 5)'
# ---------------------------------------------------------------------------

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG="$ROOT/results/validation_F1F2_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$ROOT/results"
echo "Logging to $LOG"

run_driver () {
  local method=$1 scenario=$2 projection=$3 kernel=${4:-h_w} m=${5:-1} grid=$6
  local stamp=$(date '+%H:%M:%S')
  echo "=== $stamp method=$method scenario=$scenario proj=$projection kernel=$kernel m=$m a=$grid ===" | tee -a "$LOG"
  Rscript -e "
    METHOD          <- '$method'
    SCENARIO        <- '$scenario'
    PROJECTION      <- '$projection'
    WILCOXON_KERNEL <- '$kernel'
    KERNEL          <- 'gaussian'
    M               <- ${m}L
    A_GRID          <- $grid
    N               <- ${N}L
    N_REP           <- ${N_REP}L
    NCORES          <- ${NCORES}L
    source('$ROOT/driver.R')
  " 2>&1 | tee -a "$LOG"
}

# Wilcoxon x {ito, l2, c_max_var} x m in {1, 2, 4} x {F1, F2}
for proj in ito l2 c_max_var; do
  for m in 1 2 4; do
    run_driver wilcoxon F1 "$proj" h_0 "$m" "$A_F1"
    run_driver wilcoxon F2 "$proj" h_0 "$m" "$A_F2"
  done
done

# L_p x {ito, l2, c_max_var} x p=2 x {F1, F2}
for proj in ito l2 c_max_var; do
  run_driver lp F1 "$proj" h_w 1 "$A_F1"
  run_driver lp F2 "$proj" h_w 1 "$A_F2"
done

echo "=== $(date '+%H:%M:%S') Wilcoxon (18) + L_p (6) cells complete ===" | tee -a "$LOG"

# RKHS: inline the calibrated-pipeline logic from
# run_pareto_curves_rkhs_calibrated.R, restricted to F1 and F2 so F3's
# existing results stay untouched.
echo "=== $(date '+%H:%M:%S') RKHS validation (F1, F2 only) ===" | tee -a "$LOG"
Rscript -e "
N          <- ${N}L
N_REP      <- ${N_REP}L
NCORES     <- ${NCORES}L
ALPHA      <- 0.05
SEED       <- 42L
GAMMA_GRID <- c(0, 0.25, 0.45)
A_GRID_BY_SCENARIO <- list(
  F1 = $A_F1,
  F2 = $A_F2
)
EC_GEOMETRY <- (function() {
  l_ints <- c(floor(log(N)) + 1L, N %/% 5L)
  lapply(l_ints, function(l) {
    k_star <- (N %/% 2L) - (l %/% 2L)
    c(k_star / N, (k_star + l) / N)
  })
})()
root <- '$ROOT'
res_dir <- file.path(root, 'results')
source(file.path(root, 'shared', 'integral_image.R'))
source(file.path(root, 'shared', 'weights.R'))
source(file.path(root, 'shared', 'kernels.R'))
source(file.path(root, 'shared', 'null_cache.R'))
source(file.path(root, 'shared', 'projections.R'))
source(file.path(root, 'methods', 'rkhs', 'run.R'))
suppressPackageStartupMessages(library(parallel))
T1_uniform <- function(x) rank(x, ties.method = 'average') / length(x)
T2_normal  <- function(x) qnorm(pmin(pmax(T1_uniform(x), 1e-9), 1 - 1e-9))
CELLS <- list(
  list(projection = 'ito',       kernel = 'gaussian', baseline = '',      transform = T2_normal),
  list(projection = 'l2',        kernel = 'min',      baseline = 'runif', transform = T1_uniform),
  list(projection = 'c_max_var', kernel = 'gaussian', baseline = '',      transform = T2_normal)
)
run_cell <- function(scenario, cell, a_grid) {
  source(file.path(root, 'scenarios', scenario, 'data.R'))
  generate <- get(paste0('generate_', scenario))
  proj_fn <- projection_fn(cell\$projection)
  cat(sprintf('[rkhs/%s proj=%s kernel=%s baseline=%s a_grid=(%s)] starting\n',
              scenario, cell\$projection, cell\$kernel,
              if (nzchar(cell\$baseline)) cell\$baseline else '(N(0,1))',
              paste(a_grid, collapse = ', ')))
  rows <- list()
  for (gi in seq_along(GAMMA_GRID)) {
    gamma <- GAMMA_GRID[gi]
    null_draws <- rkhs_null_distribution(
      N, gamma, kernel_name = cell\$kernel, baseline = cell\$baseline)
    cv <- critical_value(null_draws, ALPHA, m.d = 1L)
    for (eci in seq_along(EC_GEOMETRY)) {
      ec <- EC_GEOMETRY[[eci]]
      ec_length <- ec[2] - ec[1]
      for (a in a_grid) {
        set.seed(SEED + 1000L * gi + 100L * eci + as.integer(round(a * 4)))
        powers <- mclapply(seq_len(N_REP), function(rep) {
          X <- generate(N, a = a, ec.interval = ec)
          x_proj <- proj_fn(X)
          x_in <- cell\$transform(x_proj)
          stat <- rkhs_statistic(x_in, gamma, kernel = kernel_fn(cell\$kernel))
          stat > cv
        }, mc.cores = NCORES)
        power <- mean(unlist(powers))
        rows[[length(rows) + 1L]] <- data.frame(
          method = 'rkhs', scenario = scenario,
          tag = paste0(cell\$kernel,
                       if (nzchar(cell\$baseline)) paste0('_', cell\$baseline) else '',
                       '_', cell\$projection),
          n = N, m = 1L, projection = cell\$projection,
          gamma = gamma, alpha = ALPHA, N_rep = N_REP,
          ec_start = ec[1], ec_end = ec[2], ec_length = ec_length,
          a = a, power = power, critical_value = cv)
      }
    }
  }
  out <- do.call(rbind, rows)
  out_csv <- file.path(res_dir, sprintf('rkhs_%s_%s.csv', scenario, out\$tag[1]))
  write.csv(out, out_csv, row.names = FALSE)
  cat(sprintf('  wrote %s (%d rows)\n', out_csv, nrow(out)))
}
for (scenario in names(A_GRID_BY_SCENARIO)) {
  a_grid <- A_GRID_BY_SCENARIO[[scenario]]
  for (cell in CELLS) run_cell(scenario, cell, a_grid)
}
cat(sprintf('\n=== %s : 6 RKHS cells (2 scenarios x 3 proj/kernel) complete ===\n',
            format(Sys.time(), '%H:%M:%S')))
" 2>&1 | tee -a "$LOG"

echo "=== $(date '+%H:%M:%S') Validation re-run complete. Aggregate next: ==="
echo "  Rscript $ROOT/aggregate_pareto_curves.R"
