# fdachange

[![R-CMD-check](https://github.com/karbar359/fdachange/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/karbar359/fdachange/actions/workflows/R-CMD-check.yaml)

Epidemic change-point tests for functional data in R.

An *epidemic* (changed-segment) alternative means the sample
`X_1, ..., X_n` is i.i.d. except on an unknown window
`I* = {k*+1, ..., k*+l*}` where the distribution changes and then changes
back. `fdachange` implements four weighted scan tests for this alternative,
with Hoelder-type weighting that tunes sensitivity toward short windows,
ready-made critical values, p-values, and change-window localisation:

| Test | Function | Sensitive to | Reference |
|---|---|---|---|
| Wilcoxon-type antisymmetric-kernel scan | `wilcoxon_statistic()` | location (`h_w`) or scale/tail (`h_0`) changes | Bartkus and Rackauskas (2025), thesis ch. 3 |
| Lp empirical-process scan | `lp_statistic()` | general distributional change | Bartkus and Rackauskas (2026), thesis ch. 4 |
| RKHS kernel mean-embedding scan | `rkhs_statistic()` | general distributional change | Bartkus and Rackauskas (2026), thesis ch. 5 |
| Birbilas CUSUM-type functional statistic | `far1_statistic()` | functional mean jumps under FAR(1) dependence | Birbilas (2025), thesis ch. 6 |

Functional observations enter either directly (`far1_statistic`) or
through a functional-to-scalar projection (`project_l2()`, `project_ito()`,
`project_c_max_var()`, ...), optionally preceded by curve reconstruction
(`reconstruct()`). The dissertation's benchmark scenarios ship as
`generate_s1()` ... `generate_f5()` (F5: Wiener random-walk curves with a
spike-type change, restored from the paper-1 draft). A shipped dataset,
`telco_far1` (eleven simulated FAR(1) customer clusters, see
`?telco_far1`), serves as the worked application.

## Installation

```r
# install.packages("remotes")
remotes::install_github("karbar359/fdachange")
```

While the repository is private, set `GITHUB_PAT` to a token with repo
scope first (e.g. `Sys.setenv(GITHUB_PAT = gitcreds::gitcreds_get()$password)`
or a fine-grained PAT); once public, the plain call works.

The core tests run on base R alone. The basis-expansion projections
(`c_max_var`, `wiener`, `pc1`, `pc_d`) and the least-squares/penalised
reconstructions additionally need `fda` and `fda.usc`:

```r
install.packages(c("fda", "fda.usc"))
```

## Quick start

```r
library(fdachange)

## scalar series with an epidemic mean shift
set.seed(1)
x <- rnorm(200)
x[81:140] <- x[81:140] + 1

epidemic_test(x, method = "wilcoxon", gamma = 0.25)
#>
#>     Epidemic change scan test (Wilcoxon-type antisymmetric kernel)
#>
#> data:  x
#> T = 4.1877, gamma = 0.25, p-value < 2.2e-16
#> alternative hypothesis: an epidemic (changed then restored) segment
#>
#> critical value (alpha = 0.05): 2.5379  [30000 null draws]
#> detected window: observations 81..133 (fractions 0.400..0.665)

## functional data: project, then test (the ecdf-based Lp path is
## calibration-safe for any projected marginal)
X <- generate_f1(100, a = 3)          # curves with heavy-tailed change window
epidemic_test(project_l2(X), method = "lp", gamma = 0.25)

## functional data without projection (Birbilas's statistic)
Xf <- generate_f4(100, a = 2)         # FAR(1) functional time series
epidemic_test(Xf, method = "far1", gamma = 0.4)
```

## Walkthroughs

| Walkthrough | What it shows |
|---|---|
| [`vignettes/fdachange.Rmd`](vignettes/fdachange.Rmd) | The R tour: the benchmark scenarios, every projection on the same curves, calibration, a miniature power study and a reconstruct-project-test application. Read it in R as `vignette("fdachange")` after installing with `build_vignettes = TRUE`. |
| [`python/examples/fdachange.ipynb`](python/examples/fdachange.ipynb) | The same tour in Python, executed, ending with the `telco_far1` application. |
| [`reproduction/notebooks/simulations.md`](reproduction/notebooks/simulations.md) | The thesis power study, the restored F5 spike study, the Lp exponent, asymptotic vs finite-sample calibration, and what the statistic audit changed. |
| [`reproduction/notebooks/litgrid.md`](reproduction/notebooks/litgrid.md) | The LITGRID application: weekly Lithuanian balancing-price curves, legacy vs fixed pipeline, the corrected result tables. |
| [`reproduction/notebooks/telco.md`](reproduction/notebooks/telco.md) | The telco application: all eleven `telco_far1` clusters scanned, with exact calibration under FAR(1) dependence. |

The three `reproduction` notebooks render from committed result tables,
so their figures never depend on rerunning the multi-hour sweeps.

## Critical values and null distributions

Every test needs the null distribution of its scan statistic. The package
resolves it automatically, in this order:

1. **Closed form** where it exists: the Wilcoxon scan at `gamma = 0` is
   calibrated by Kuiper's law (`kuiper_cdf()`, `kuiper_quantile()`),
   exactly.
2. **Shipped Monte-Carlo caches** (`inst/extdata`, about 2 MB): the
   dissertation's null draws for all published (gamma, kernel, p) cells;
   accessed via `wilcoxon_null_draws()`, `lp_null_draws()`,
   `rkhs_null_draws()`, `far1_null_draws()`.
3. **Your own draws**: for uncached cells, run the matching simulator
   (`simulate_wilcoxon_null()`, `simulate_lp_null()`,
   `simulate_rkhs_null()`, `simulate_far1_null()`) and either store
   the result with `save_null_draws()` (picked up automatically from then
   on) or pass it directly as `epidemic_test(..., null = draws)`.

Critical values are type-7 empirical quantiles of these draws
(`critical_value()`), with an exact adjustment for split-sample
aggregation (`scan_over_splits()`, `splits =` argument).

One calibration caveat: the Wilcoxon scan and the ecdf-path Lp scan are
rank-based, so their cached nulls apply to any continuous data. The RKHS
caches are tied to the H0 sampling distribution named by their `baseline`
tag; when your (projected) data has a different null marginal, simulate
the null under your H0 (`simulate_rkhs_null(rdist = ...)`) instead of
relying on a mismatched cache.

## Provenance and verification

The implementations follow the statistics exactly as defined in the
dissertation. They were extracted from the thesis codebase after a full
statistic-level audit (September 2026) that brute-force verified every
formula; the pre-audit research code contained implementation defects
(most notably an inverted weight direction in the Wilcoxon scan) whose
printed consequences are documented in the dissertation's erratum
material. This package ships only the verified, corrected family.

The test suite enforces that permanently:

- brute-force implementations of each printed formula agree with the
  package functions to `1e-10`;
- 276 golden values captured from the audited original code (statistics,
  scan windows, critical values from every shipped cache cell, RNG
  streams of every scenario) agree to `1e-12`;
- the shipped null caches are converted from the original Monte-Carlo
  CSVs with element-wise `identical()` assertions, so critical values are
  bit-identical to the original pipeline.

The rerun pipeline behind those checks, including the vendored
pre-audit code, lives under [`reproduction/`](reproduction/): its
[`README.md`](reproduction/README.md) describes the layout and the
verification gates, [`SEEDS.md`](reproduction/SEEDS.md) the RNG
discipline that makes every run core-count independent.

## Python implementation

A Python twin lives under [`python/`](python/) with a one-to-one API
(same snake_case names), numpy/scipy internals, and null caches converted
from the same frozen sources, so critical values are bit-identical. The R
package is the reference implementation; the two are held equal by a
cross-language parity suite (statistics on identical inputs at 1e-10,
critical values at 1e-12, reconstruction outputs against R fixtures).

```sh
pip install "git+https://github.com/karbar359/fdachange#subdirectory=python"
```

`skfda.FDataGrid` inputs are accepted everywhere via the optional
scikit-fda extra (`fdachange[fda]`), which also provides the penalized
reconstruction arm; see [`python/README.md`](python/README.md) for the
API table and the ecosystem notes (scikit-fda, ruptures).

## Citing

If you use this package, please cite the dissertation and the
corresponding article:

- Bartkus, K. (2026). *Changed Segment Tests on Functional Data*.
  PhD dissertation, Vilnius University. doi:10.15388/vu.thesis.950.
- Bartkus, K., Rackauskas, A. (2025). Changed segment tests via norm and
  Brownian projection. *Lithuanian Mathematical Journal*, 65(4), 481-499.
  doi:10.1007/s10986-025-09699-7. (Wilcoxon-type scan.)
- Bartkus, K., Rackauskas, A. (2026). Nonparametric changed segment
  detection in functional data. *Nonlinear Analysis: Modelling and
  Control*, 31(1), 194-211. doi:10.15388/namc.2026.31.44489. (Lp and RKHS
  scans.)
- Birbilas, A. (2025). *Functional Data Analysis in Telecommunications*.
  PhD dissertation, Vilnius University. doi:10.15388/vu.thesis.841.
  (For `far1_statistic()`.)

See `citation("fdachange")`.

## Contributing

Issues and pull requests are welcome; the checks to run are listed in
[`CONTRIBUTING.md`](CONTRIBUTING.md). Release notes: [`NEWS.md`](NEWS.md).

## License

MIT. See [`LICENSE.md`](LICENSE.md).
