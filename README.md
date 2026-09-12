# fdachange

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
| Birbilas CUSUM-type functional statistic | `birbilas_statistic()` | functional mean jumps under FAR(1) dependence | Birbilas (2025), thesis ch. 6 |

Functional observations enter either directly (`birbilas_statistic`) or
through a functional-to-scalar projection (`project_l2()`, `project_ito()`,
`project_c_max_var()`, ...), optionally preceded by curve reconstruction
(`reconstruct()`). The dissertation's benchmark scenarios ship as
`generate_s1()` ... `generate_f4()`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("karbar359/fdachange")
```

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
#>  Epidemic change scan test (Wilcoxon-type antisymmetric kernel)
#>  T = ..., p-value = ...
#>  critical value (alpha = 0.05): ...
#>  detected window: observations ...

## functional data: project, then test (the ecdf-based Lp path is
## calibration-safe for any projected marginal)
X <- generate_f1(100, a = 3)          # curves with heavy-tailed change window
epidemic_test(project_l2(X), method = "lp", gamma = 0.25)

## functional data without projection (Birbilas's statistic)
Xf <- generate_f4(100, a = 2)         # FAR(1) functional time series
epidemic_test(Xf, method = "birbilas", gamma = 0.4)
```

## Critical values and null distributions

Every test needs the null distribution of its scan statistic. The package
resolves it automatically, in this order:

1. **Closed form** where it exists: the Wilcoxon scan at `gamma = 0` is
   calibrated by Kuiper's law (`kuiper_cdf()`, `kuiper_quantile()`),
   exactly.
2. **Shipped Monte-Carlo caches** (`inst/extdata`, about 2 MB): the
   dissertation's null draws for all published (gamma, kernel, p) cells;
   accessed via `wilcoxon_null_draws()`, `lp_null_draws()`,
   `rkhs_null_draws()`, `birbilas_null_draws()`.
3. **Your own draws**: for uncached cells, run the matching simulator
   (`simulate_wilcoxon_null()`, `simulate_lp_null()`,
   `simulate_rkhs_null()`, `simulate_birbilas_null()`) and either store
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

## Citing

If you use this package, please cite the dissertation and the
corresponding article:

- Bartkus, K. (2026). *Analysis of Epidemic Changes in Functional Data*.
  PhD dissertation, Vilnius University.
- Bartkus, K., Rackauskas, A. (2025). Lithuanian Mathematical Journal.
  (Wilcoxon-type scan.)
- Bartkus, K., Rackauskas, A. (2026). Nonlinear Analysis: Modelling and
  Control. (Lp and RKHS scans.)
- Birbilas, A. (2025). *Functional Data Analysis in Telecommunications*.
  PhD dissertation, Vilnius University. doi:10.15388/vu.thesis.841.
  (For `birbilas_statistic()`.)

See `citation("fdachange")`.

## License

MIT. See `LICENSE.md`.
