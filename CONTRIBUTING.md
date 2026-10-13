# Contributing

Issues and pull requests are welcome.

- Open PRs against `main`; one logical change per PR.
- Run `Rscript -e 'roxygen2::roxygenise()'` after changing docs and make
  sure `R CMD check` and `testthat::test_local()` pass (the golden-value
  suite must stay green: it pins the statistics to the audited thesis
  implementation).
- Slow distributional tests run with `FDACHANGE_RUN_SLOW=1`.
- Python changes: `pip install './python[test]'`, then
  `pytest python/tests -q`. The example notebook must still execute end
  to end: `jupyter execute python/examples/fdachange.ipynb` (needs
  matplotlib, nbclient and ipykernel).
- Changes to committed result data or to the notebook sources under
  `reproduction/notebooks/`: install the package (`R CMD INSTALL .`),
  re-render with `Rscript reproduction/notebooks/render.R` and commit
  the regenerated `.md` files and figures.

Each sub-tree documents its own layout: [`python/README.md`](python/README.md)
and [`reproduction/README.md`](reproduction/README.md).
