#!/usr/bin/env Rscript
# Aggregate the per-cell F4 driver outputs into:
#   - `results/F4_combined.csv`     -- long format, all rows from every
#                                       <method>_F4*.csv concatenated.
#   - `results/F4_summary_wide.csv` -- wide pivot keyed on (method, cell,
#                                       gamma, ec_length, n_rep) with
#                                       one column per `a` value.
# Plus a printed wide table at the canonical ec_length (l/n = 0.2) for at-
# a-glance reading.
#
# Usage:
#   Rscript reproduction/aggregate_F4.R

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- sub(key, "", args[grepl(key, args)])
  if (length(hit)) return(normalizePath(hit[1]))
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) return(normalizePath(ofile))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))
res_dir <- file.path(root, "results")

files <- list.files(res_dir,
                    pattern = "^(lp|wilcoxon|rkhs|birbilas)_F4.*\\.csv$",
                    full.names = TRUE)
files <- files[!grepl("_exact|_corrected", files)]   # corrected-statistic re-runs; see compare_exact_vs_legacy.R
if (!length(files)) stop("No F4 result CSVs in ", res_dir)
cat(sprintf("Aggregating %d CSV(s) from %s\n", length(files), res_dir))

# Long format: just rbind everything.
long <- do.call(rbind, lapply(files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$source_file <- basename(f)
  d
}))

# Build a `cell` label per row: `<method>:<tag>`. The driver's out_tag
# already encodes projection / kernel / m-split / etc., so the cell
# label inherits all of that without duplication.
long$cell <- paste(long$method, long$tag, sep = ":")

long_path <- file.path(res_dir, "F4_combined.csv")
write.csv(long, long_path, row.names = FALSE)
cat("wrote ", long_path, " (", nrow(long), " rows)\n", sep = "")

# Wide pivot: one row per (method, cell, gamma, ec_length, N_rep), one
# column per `a` value. Reshape via base R.
key_cols  <- c("method", "cell", "n", "m", "projection", "gamma",
               "ec_start", "ec_end", "ec_length", "N_rep", "alpha")
wide <- reshape(long[, c(key_cols, "a", "power")],
                idvar = key_cols, timevar = "a", direction = "wide",
                v.names = "power")
# rename "power.X" -> "a=X"
a_cols <- grep("^power\\.", names(wide), value = TRUE)
names(wide)[match(a_cols, names(wide))] <- sub("^power\\.", "a=", a_cols)

# Sort for readable output.
wide <- wide[order(wide$method, wide$cell, wide$ec_length, wide$gamma), ]

wide_path <- file.path(res_dir, "F4_summary_wide.csv")
write.csv(wide, wide_path, row.names = FALSE)
cat("wrote ", wide_path, " (", nrow(wide), " rows)\n", sep = "")

# Print the canonical-ec-length slice (l/n = 0.2) as a quick-look table.
canon <- wide[abs(wide$ec_length - 0.2) < 1e-3, , drop = FALSE]
if (nrow(canon)) {
  cat("\n--- Headline table at ec_length = 0.2 (l = n/5 = 40 at n=200) ---\n")
  display_cols <- c("cell", "gamma",
                    grep("^a=", names(canon), value = TRUE))
  options(width = 220)
  print(canon[, display_cols], row.names = FALSE)
}
