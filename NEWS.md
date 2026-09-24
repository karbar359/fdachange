# fdachange (development version)

- Breaking: the functional FAR(1)-type statistic's API is renamed from the
  author name to `far1`: `birbilas_statistic()` -> `far1_statistic()`,
  `birbilas_scan()` -> `far1_scan()`, `birbilas_betas` -> `far1_betas`,
  `birbilas_null_draws()` -> `far1_null_draws()`,
  `birbilas_critical_value()` -> `far1_critical_value()`,
  `simulate_birbilas_null()` -> `simulate_far1_null()`, and the method
  string `"birbilas"` -> `"far1"` in `epidemic_test()` and
  `save_null_draws()`. The shipped null cache is now
  `inst/extdata/far1_null.rds`; a user cache saved as `birbilas_null.rds`
  in `null_cache_dir()` must be renamed to `far1_null.rds`.

# fdachange 0.1.2

- New vignette "Detecting epidemic changes in functional data"
  (`vignette("fdachange")`): quick start, the benchmark scenarios
  including F5's spike function, every projection visualized on the same
  curves, the calibration-resolution story, a miniature power study with
  the full committed study's verdict, and an application-shaped
  reconstruct-project-test walkthrough with the corrected LITGRID
  results quoted.

# fdachange 0.1.1

- New scenario `generate_f5()`: i.i.d. Wiener (random-walk) curves with a
  spike `a * f_alpha(t)` added on the epidemic window, where `f_alpha` is
  the dissertation's spike function
  `(1 - |t - 1/2|^alpha)^(1/alpha)` (`spike_alpha` argument, default
  0.45 as in the Faber-Schauder smoothing figure). Restored from the
  paper-1 draft's scenario (S2) (`p1.TD.S2.spike` / `wiener_TD2_spike` in
  the archived research code, which used a tent-power bump), cut before
  publication and lost in the original extraction; it is the experiment
  behind the thesis's claim that point-evaluation projections suit
  spike-type changes. `scenario_cdf("F5")` returns the ecdf transform
  like the other functional scenarios.

# fdachange 0.1.0

Initial release: the epidemic change-point test family of Bartkus (2026),
extracted from the thesis codebase after the September 2026 statistic
audit.

- Four scan statistics with localisation: `wilcoxon_statistic()` /
  `wilcoxon_scan()`, `lp_statistic()` / `lp_scan()`, `rkhs_statistic()` /
  `rkhs_scan()`, `birbilas_statistic()` / `birbilas_scan()`.
- `epidemic_test()`: htest-style wrapper with p-values, critical values
  and the detected window; split-sample aggregation via `splits =`.
- Null calibration: closed-form Kuiper law at `gamma = 0` (Wilcoxon),
  shipped Monte-Carlo caches for all published cells, on-demand
  simulators plus a user-side cache (`save_null_draws()`,
  `null_cache_dir()`) for everything else.
- Functional-data support: projections (`project_l2()`, `project_ito()`,
  `project_c_max_var()`, `project_wiener()`, `project_pc1()`,
  `project_pc_d()`), curve reconstruction (`reconstruct()`: Faber-Schauder
  interpolation, least squares, penalised B-splines), and the benchmark
  scenarios `generate_s1()` ... `generate_f4()`.

Lineage note: this package ships the corrected statistic family. The
implementations equal brute-forced versions of the dissertation's printed
formulas to machine precision and are pinned to the audited original code
by a golden-value test suite. The pre-audit research implementations that
produced some published tables are not part of the package; they are kept
in the repository's `reproduction/` tree solely to make those runs
reproducible.
