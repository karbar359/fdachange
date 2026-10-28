# RNG discipline and seeds

Every stochastic entry point sets `RNGkind("L'Ecuyer-CMRG")` and a fixed
seed, and parallelises either with `parallel::mclapply(mc.set.seed = TRUE)`
(driver, smoke) or with one explicit `nextRNGStream` substream per draw
(the null-cache generators), so each fork receives its own deterministic
RNG substream.

## Per-script seeds

| Script | Seed | Stream structure |
|---|---|---|
| `driver.R` | `SEED = 42` | one L'Ecuyer substream per (gamma, ec, a) task, assigned in `expand.grid(gamma, ec_idx, a)` order via `mc.preschedule = FALSE` |
| `smoke.R` | `set.seed(42)` | reset inside each (cell, a) evaluation, so cells are order-independent |
| `stability.R` | `SEED = 42` | `set.seed(SEED + cell_idx)` per fixed-sample cell: fully order-independent |
| `run_pareto_curves_rkhs_*.R`, `run_thesis_power_F1F2_rkhs_*.R` | `SEED = 42` | `set.seed(SEED + 1000*gamma_idx + 100*ec_idx + round(4*a))` per cell: order-independent |
| `verify_statistics.R` | `SEED = 7` | serial |
| `applications/litgrid/*.R` | `SEED = 3` | serial; the 10 Ito draws are drawn once and shared across methods and gammas |
| `generators/*.R` (null caches) | `SEED = 2` | one L'Ecuyer substream per draw, chained with `nextRNGStream` in `mclapply_with_progress`; the parent moves to the stream after the last draw, so cells are order-dependent but core-count independent |
| package golden capture (`data-raw/capture_golden_values.R`) | `seed 7 + n` per data block | serial |

## Core-count independence

`driver.R` results are independent of `NCORES`: with
`mc.preschedule = FALSE` jobs are forked in task order and each fork takes
the next L'Ecuyer substream, so the stream-to-task assignment never depends
on how many cores run. Verified empirically on 2026-09-12: the
`wilcoxon/S1 n=400 corrected` sweep at `NCORES = 3` and `NCORES = 20`
produced byte-identical CSVs, both equal to the original pipeline's output.

The null-cache generators are core-count independent as well:
`mclapply_with_progress` assigns draw `i` the `i`-th `nextRNGStream`
substream, whatever `ncores` or `batch_size`. Until 2026-10 it relied on
`mc.set.seed = TRUE` plus a `runif(1)` between batches, which re-forked
the same substreams every second batch shifted by one normal draw: the
patched `gamma = 0` Wilcoxon column (acf 0.994 at lag 160), the RKHS
`SD_h_n2000_corrected` and far1 `E_q_gamma` caches (acf 0.84 to 0.92) and
the finite F4 caches `E_n_beta_gamma_finite*` and `SD_h_F4_finite*` (acf
0.66 to 0.84) were produced that way and held about 100 to 160 distinct
draws each. All of them were regenerated with the fixed helper on
2026-10-28 (`data-raw/cache-src/`; the gamma = 0 rows carry the CSV tag
`port_wilcoxon_null_N10000_MC30000_independent`), and the committed tables
in `data/` that read them were refreshed with
`sweeps/refresh_regenerated_nulls.sh`. The Wilcoxon `gamma > 0` and Lp
caches are the original pipeline's sequential runs and the legacy RKHS
cache `SD_h_n2000_dist_t` is the thesis's own 300-draw artefact; none of
them was affected. Regenerated caches match the frozen ones
distributionally (see `generators/validate.R`), not draw for draw.

## Corrected-family sweep configuration

The corrected re-run of the simulation chapter
(`sweeps/run_corrected_n400.sh`) used: `EC_EXACT = TRUE` (exact-length
change windows), `N = 400`, `SEED = 42`, `N_REP = 100` except F1/F2 cells
at the printed `N_REP = 20`, `NCORES = 20` (immaterial, see above), and the
default grids `GAMMA_GRID = c(0, 0.25, 0.45)`,
`A_GRID = c(0, 0.25, 0.5, 1, 2, 4, 8)` for scalar scenarios with
scenario-specific a-grids for F1-F3.

## Environment

R 4.3.3 on x86_64 Linux; package versions pinned in `renv.lock`. The
statistics themselves are pure R arithmetic (no BLAS-order-sensitive
reductions), so results are expected to be bit-stable across machines with
identical package versions; the golden-value suite in the package tests
enforces this.
