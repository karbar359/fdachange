#!/usr/bin/env Rscript
# Aggregate per-cell driver outputs for the Pareto-curve scenarios
# (F1, F2, F3) into a long CSV (`pareto_curves_combined.csv`) and a
# wide CSV (`pareto_curves_summary_wide.csv`) mirroring the F4
# aggregator output format.

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
                    pattern = "^(lp|wilcoxon|rkhs)_F[123].*\\.csv$",
                    full.names = TRUE)
files <- files[!grepl("_exact|_corrected", files)]   # corrected-statistic re-runs; see compare_exact_vs_legacy.R
if (!length(files)) stop("No F1/F2/F3 result CSVs in ", res_dir)
cat(sprintf("Aggregating %d CSV(s) from %s\n", length(files), res_dir))

long <- do.call(rbind, lapply(files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$source_file <- basename(f)
  d
}))
long$cell <- paste(long$method, long$tag, sep = ":")

long_path <- file.path(res_dir, "pareto_curves_combined.csv")
write.csv(long, long_path, row.names = FALSE)
cat("wrote ", long_path, " (", nrow(long), " rows)\n", sep = "")

key_cols  <- c("method", "scenario", "cell", "n", "m", "projection", "gamma",
               "ec_start", "ec_end", "ec_length", "N_rep", "alpha")
wide <- reshape(long[, c(key_cols, "a", "power")],
                idvar = key_cols, timevar = "a", direction = "wide",
                v.names = "power")
a_cols <- grep("^power\\.", names(wide), value = TRUE)
names(wide)[match(a_cols, names(wide))] <- sub("^power\\.", "a=", a_cols)
wide <- wide[order(wide$scenario, wide$method, wide$cell,
                   wide$ec_length, wide$gamma), ]
wide_path <- file.path(res_dir, "pareto_curves_summary_wide.csv")
write.csv(wide, wide_path, row.names = FALSE)
cat("wrote ", wide_path, " (", nrow(wide), " rows)\n", sep = "")
