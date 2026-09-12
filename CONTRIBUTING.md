# Contributing

Issues and pull requests are welcome.

- Open PRs against `main`; one logical change per PR.
- Run `Rscript -e 'roxygen2::roxygenise()'` after changing docs and make
  sure `R CMD check` and `testthat::test_local()` pass (the golden-value
  suite must stay green: it pins the statistics to the audited thesis
  implementation).
- Slow distributional tests run with `FDACHANGE_RUN_SLOW=1`.
