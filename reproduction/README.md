# Reproduction of the dissertation's computational results

This directory reruns the simulation and application pipeline of Bartkus
(2026), *Changed Segment Tests on Functional Data* (PhD dissertation,
Vilnius University, doi:10.15388/vu.thesis.950), on top of the in-tree
`fdachange` package. It exists so that
every number the pipeline produced - including the pre-audit ("legacy")
runs behind some published tables - can be regenerated and compared.

It is excluded from the built R package (`.Rbuildignore`).

## Layout

| Path | What it is |
|---|---|
| `R/legacy/` | The original implementation files, vendored verbatim: all four statistic families (legacy, exact AND corrected variants), weights, kernels, null-cache loaders, projections, reconstruction. The pre-audit variants exist only here; the package ships only the corrected family. |
| `R/compat.R` | Original-pipeline conventions (`generate_S1`, `cdf_F2`, `EC_EXACT`, ...) mapped onto the package, whose generators are bit-identical to the originals. |
| `R/legacy_env.R` | `legacy_env()` loads the vendored files into an isolated environment; `legacy_attach()` into the global one (for pure-legacy scripts). |
| `driver.R` | The simulation driver: one (method, scenario) sweep per run; variant flags pick `fdachange::` (corrected) or `legacy$` (pre-audit). |
| `verify_statistics.R` | The audit gate: brute-forced printed formulas vs both families, Kuiper law vs printed Table 3.1, scenario-level checks, package == vendored-corrected identity. Run this first. |
| `smoke.R`, `stability.R`, `mc_calibration_demo.R`, `benchmark_scaling.R` | End-to-end smoke table, projection-stability study, Monte-Carlo calibration demo, timing/memory limits. |
| `compare_exact_vs_legacy.R` | Pairs re-run family CSVs with legacy CSVs and writes the committed `data/{exact,corrected}_vs_legacy.csv`. |
| `sweeps/` | Shell wrappers that batch `driver.R` over the thesis grids; `run_corrected_n400.sh` is the corrected-family re-run of the whole simulation chapter. |
| `generators/` | Monte-Carlo null-cache generators + `validate.R` (distributional check against the frozen caches). |
| `aggregators/` | Table/figure builders; LaTeX fragments and PNGs go to `output/` (gitignored). |
| `notebooks/` | Rendered walkthroughs (view the `.md` files on GitHub): `simulations` visualizes the committed power/calibration/audit tables, `litgrid` and `telco` the two applications. Rebuild with `Rscript notebooks/render.R`; rendered in CI against the committed data. Python users: package usage and the telco application live in `python/examples/fdachange.ipynb`. |
| `applications/litgrid/` | The LITGRID balancing-price application, including `reproduce_tables.R` (asserts the printed Tables 7.1/7.2 within print tolerance) and `download_litgrid.R`. |
| `data/` | Committed tidy outputs (the canonical result CSVs the thesis figures used). |
| `results/`, `output/` | Gitignored working state. |

## Inputs

- **Null caches**: `../data-raw/cache-src/` (committed; the frozen originals).
  Override with `REPRO_CACHE_DIR`.
- **LITGRID raw data** (~119 MB, gitignored): fetch with
  `Rscript applications/litgrid/download_litgrid.R` (GitHub release asset
  `litgrid-raw-v1`, checksum-verified against the committed `SHA256SUMS`).
  Override with `LITGRID_RAW_DIR`.
- **Thesis tex fragments** (only for `aggregators/assemble_thesis_power.R`'s
  printed-value extraction): point `THESIS_TEX_DIR` at the thesis source's
  simulation-chapter fragments; without them, the committed
  `data/thesis_power_values.csv` remains the canonical artifact.

## The restored F5 scenario and its power study

`fdachange::generate_f5()` restores the paper-1 draft's scenario (S2)
(archive ids `p1.TD.S2.spike` / `wiener_TD2_spike`, cut before
publication): i.i.d. Wiener curves with a spike `a f_alpha(t)` added on
the epidemic window, where `f_alpha(t) = (1 - |t - 1/2|^alpha)^(1/alpha)`
is the dissertation's spike function (default `alpha = 0.45`, as in the
Faber-Schauder smoothing figure; the archived draft/code used a
tent-power bump instead). It appears in no printed thesis table, but it
is the experiment behind the thesis's claim that point-evaluation
projections suit spike-type changes
(`ch2_dimensionality/sec_point_eval.tex`).

`sweeps/run_f5_spike_power.sh` runs that experiment (corrected Wilcoxon
h_0 and Lp p = 2 over the projections c_max_var / l2 / ito, n = 400,
N_rep = 100, seed 42); `aggregators/aggregate_f5_spike.R` combines the
outputs into the committed `data/f5_spike_power.csv` and prints the
projection comparison.

**Result (2026-09-16 run):** the point-evaluation claim holds in the
weak-signal regime but not uniformly. For the Lp scan, c_max_var clearly
leads where the spike is small (long window, a = 0.5-1: power 0.16-0.57
vs 0.06-0.28 for l2) and has the best mean power over all alternatives
(0.46 vs 0.41 ito, 0.40 l2). At large amplitudes the norm projection
catches up (l2 reaches 1.00 where c_max_var sits at 0.84-0.99) and
c_max_var is visibly non-monotone in a - its single randomly chosen
variance-critical point becomes unstable when the spike dominates the
variance function. For the Wilcoxon h_0 scan the ordering reverses on
average (l2 0.34, c_max_var 0.30, ito 0.24). Empirical size stays near
the nominal 0.05 in almost all cells (max 0.16 in one lp/ito cell,
within a few MC standard errors at N_rep = 100). In short: point
evaluations are the most informative projection for *faint* spike-type
changes, exactly the regime the thesis argument concerns, but they are
not uniformly dominant.

## Verification gates (all green as of 2026-09-12)

1. `Rscript verify_statistics.R` - every corrected statistic equals its
   brute-forced printed formula at 1e-10, every legacy defect mechanism
   reproduces, and the package equals the vendored corrected variants bit
   for bit.
2. Driver bit-reproduction - the corrected-family gate cell
   (`wilcoxon/S1, n = 400, seed 42, EC_EXACT`) regenerated through this
   driver is byte-identical to the original pipeline's
   `wilcoxon_S1_n400_h_w_beta0_corrected.csv`, at any core count (see
   `SEEDS.md`).
3. `Rscript applications/litgrid/reproduce_tables.R` - the printed
   dissertation Tables 7.1 and 7.2 reproduce within print tolerance
   (exits non-zero otherwise).

## Typical runs

```sh
Rscript verify_statistics.R                       # < 1 min, run first
Rscript smoke.R                                   # end-to-end table, minutes
bash sweeps/run_corrected_n400.sh                 # corrected re-run, hours
Rscript -e "FAMILY <- 'corrected'; source('compare_exact_vs_legacy.R')"
Rscript applications/litgrid/reproduce_tables.R   # printed-table assertion
Rscript notebooks/render.R                        # re-render the notebooks
```

Seeds, RNG streams and the core-count-independence proof: `SEEDS.md`.

## Environment

Package versions are recorded in `renv.lock` (R 4.3.3). The scripts run
against your regular library; to recreate the recorded environment,
install the packages from the lockfile, e.g.

```r
renv::restore(lockfile = "reproduction/renv.lock")
```

The lockfile is a record, not an activated renv project: there is
deliberately no `.Rprofile` here, so the scripts keep seeing the library
that has `fdachange` installed.
