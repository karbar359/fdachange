# fdachange (Python)

Epidemic change-point tests for functional data - the Python twin of the
[fdachange R package](https://github.com/karbar359/fdachange), which is
the **reference implementation**. The two are held equal by a
cross-language parity suite: statistics, scan windows and critical values
are asserted in pytest against values computed by the R package on
identical inputs, and the packaged null caches are converted from the
same frozen Monte-Carlo sources, so critical values are bit-identical.

## Install

```sh
pip install "git+https://github.com/karbar359/fdachange#subdirectory=python"
# with scikit-fda interop and the penalized reconstruction arm:
pip install "fdachange[fda] @ git+https://github.com/karbar359/fdachange#subdirectory=python"
```

Core dependencies: numpy, scipy. The `[fda]` extra adds scikit-fda.

## Quick start

```python
import fdachange as fc

x = fc.generate_s1(200, a=1, rng=1)            # epidemic mean shift
res = fc.epidemic_test(x, "wilcoxon", gamma=0.25)
print(res)                                     # statistic, p-value, CV, window
res.reject, res.interval

X = fc.generate_f5(100, a=4, rng=2)            # Wiener curves + spike
fc.epidemic_test(fc.project_l2(X), "lp", gamma=0.25)

telco = fc.load_telco_far1()                   # shipped FAR(1) dataset
fc.epidemic_test(fc.project_l2(telco["cl2"].y.T), "wilcoxon", gamma=0.25)
```

Functional inputs are `(grid, n_observations)` arrays - one curve per
column, mirroring the R convention - or `skfda.FDataGrid` objects (any
matrix-consuming function accepts both).

For a worked tour - scenarios, projections visualized, calibration, a
mini power study and an application on the shipped dataset - see the
executed example notebook
[`examples/fdachange.ipynb`](examples/fdachange.ipynb), the counterpart
of the R package's vignette (executed in CI, so it stays current).

## API

Mirrors the R package one-to-one (all names are identical):

| Area | Functions |
|---|---|
| Statistics + localization | `wilcoxon_statistic/scan`, `lp_statistic/scan`, `rkhs_statistic/scan`, `far1_statistic/scan`, `epidemic_test` |
| Calibration | `kuiper_cdf/quantile` (exact gamma = 0 law), `*_null_draws`, `*_critical_value`, `simulate_*_null`, `save_null_draws`, `null_cache_dir` |
| Scenarios | `generate_s1..s5`, `generate_f1..f5`, `scenario_cdf` |
| Projections | `project_l2/ito/wiener/c_max_var/pc1/pc_d`, `projection_fn`, `projection_multiplicity` |
| Reconstruction | `reconstruct` (`fs_interp` = Faber-Schauder truncation proper; `fs_ls` = least squares onto the same polygonal space; `bspline_pen` via scikit-fda) |
| Aggregation | `scan_over_splits`, `make_run` |
| Data | `load_telco_far1` (the R package's `telco_far1` dataset, bit-identical) |

## Differences from the R package

- **RNG**: numpy `Generator` (pass `rng=` as a Generator, a seed, or
  `None`); draw streams are not reproducible across languages. Parity is
  established on fixed inputs, not on streams.
- **`bspline_pen`** needs the `[fda]` extra (scikit-fda `BasisSmoother`).
- Heavy null-cache regeneration is nicer in R (parallel simulators);
  the Python simulators are serial reference implementations.

## Ecosystem

[scikit-fda](https://fda.readthedocs.io) provides the general FDA
toolbox (representations, smoothing, FPCA, registration) and its
`FDataGrid` plugs directly into this package;
[ruptures](https://centre-borelli.github.io/ruptures-docs/) covers
classical change-point detection. Neither implements epidemic-window scan
statistics with Hoelder weighting and distribution-free calibration -
the niche fdachange fills.

## Citing

Cite the dissertation and articles listed by `citation("fdachange")` in
the R package (see the repository root README).
