#!/usr/bin/env Rscript
# Export data/telco_far1.rda to full-precision (%.17g) CSVs for the Python
# packer (build_telco_data.py), one set of files per cluster into the
# gitignored python/tools/telco-src/. The committed artifact is the .npz
# the packer builds; these CSVs are the intermediate hand-off.
#
# Run from anywhere:  Rscript python/tools/export_telco.R

here <- dirname(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
load(normalizePath(file.path(here, "..", "..", "data", "telco_far1.rda")))
out_dir <- normalizePath(file.path(here, "telco-src"), mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE)

g17 <- function(x) formatC(x, format = "g", digits = 17)
write_mat <- function(m, name) {
  write.table(matrix(g17(m), nrow = nrow(m)), file.path(out_dir, name),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
write_vec <- function(x, name, fmt = g17) {
  writeLines(fmt(x), file.path(out_dir, name))
}

for (nm in names(telco_far1)) {
  Y <- telco_far1[[nm]]
  write_mat(unclass(Y), sprintf("%s_y.csv", nm))
  write_vec(attr(Y, "change_months"), sprintf("%s_change_months.csv", nm),
            as.character)
  write_vec(attr(Y, "multiplier"), sprintf("%s_multiplier.csv", nm))
  write_vec(attr(Y, "jump"), sprintf("%s_jump.csv", nm))
  out <- attr(Y, "outliers")
  write.table(out, file.path(out_dir, sprintf("%s_outliers.csv", nm)),
              sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write_vec(attr(Y, "outlier_months"), sprintf("%s_outlier_months.csv", nm),
            as.character)
}
cat(sprintf("wrote %d clusters to %s\n", length(telco_far1), out_dir))
