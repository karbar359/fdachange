# Reproduction of the dissertation's computational results

This directory reruns the simulation and application pipeline of Bartkus
(2026), *Analysis of Epidemic Changes in Functional Data* (Vilnius
University), on top of the in-tree `fdachange` package. It exists so that
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
```

Seeds, RNG streams and the core-count-independence proof: `SEEDS.md`.
Package environment: `renv.lock` (restore with `renv::restore()`).
