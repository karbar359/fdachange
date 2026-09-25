# Load the vendored pre-audit implementations into an isolated environment,
# so legacy names (wilcoxon_statistic, load_rkhs_null, ...) never shadow the
# fdachange package. The environment holds the complete original files -
# legacy, exact AND corrected variants plus the null-cache loaders - exactly
# as they ran in the thesis pipeline; entry-point scripts select
# `legacy$<fn>` for pre-audit variants and `fdachange::<fn>` for the
# corrected family.
#
# `projections = TRUE` additionally sources the legacy projections and
# reconstruction files (attaches fda / fda.usc); the ported entry points use
# the package versions, so it is off by default.

.legacy_files <- function(projections) {
  files <- c("integral_image.R", "weights.R", "kernels.R", "null_cache.R",
             "progress.R", "splits.R",
             "method_wilcoxon.R", "method_lp.R", "method_rkhs.R",
             "method_far1.R")
  if (projections) files <- c(files, "projections.R", "reconstruction.R")
  files
}

legacy_env <- function(repro_root, projections = FALSE) {
  env <- new.env(parent = globalenv())
  for (f in .legacy_files(projections)) {
    source(file.path(repro_root, "R", "legacy", f), local = env)
  }
  env
}

# Source the vendored files into the GLOBAL environment, old-pipeline style.
# For the pure-legacy analysis scripts (litgrid application, generators)
# whose bodies then run unchanged. Do not mix with package calls in the same
# script: the legacy names shadow the package exports.
legacy_attach <- function(repro_root, projections = TRUE) {
  for (f in .legacy_files(projections)) {
    source(file.path(repro_root, "R", "legacy", f))
  }
  invisible(TRUE)
}
