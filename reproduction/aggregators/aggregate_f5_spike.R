#!/usr/bin/env Rscript
# Combine the F5 (Wiener + spike) power-study CSVs written by
# sweeps/run_f5_spike_power.sh into the committed
# data/f5_spike_power.csv, and print the projection comparison that the
# study exists for: does the point-evaluation projection (c_max_var) beat
# the integral projections (l2, ito) on a spike-type change?
#
# Usage:  Rscript reproduction/aggregators/aggregate_f5_spike.R

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) return(normalizePath(hit[1]))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))

files <- list.files(file.path(root, "results"),
                    pattern = "^(wilcoxon|lp)_F5_n400_.*_corrected_.*\\.csv$",
                    full.names = TRUE)
if (!length(files)) stop("No F5 result CSVs in results/; run sweeps/run_f5_spike_power.sh first.")

out <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
out <- out[order(out$method, out$projection, out$ec_length, out$gamma, out$a), ]

out_path <- file.path(root, "data", "f5_spike_power.csv")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_path, row.names = FALSE)
cat(sprintf("wrote %s (%d rows from %d files)\n\n", out_path, nrow(out), length(files)))

# ---- the comparison of record ----------------------------------------------
cat("Power by (method, window, gamma, a) across projections:\n")
alt <- out[out$a > 0, ]
wide <- reshape(alt[, c("method", "ec_length", "gamma", "a", "projection", "power")],
                idvar = c("method", "ec_length", "gamma", "a"),
                timevar = "projection", direction = "wide")
names(wide) <- sub("^power\\.", "", names(wide))
wide <- wide[order(wide$method, wide$ec_length, wide$gamma, wide$a), ]
print(format(wide, digits = 2), row.names = FALSE)

cat("\nSize (a = 0):\n")
size <- out[out$a == 0, c("method", "projection", "ec_length", "gamma", "power")]
print(format(reshape(size, idvar = c("method", "ec_length", "gamma"),
                     timevar = "projection", direction = "wide"),
             digits = 2), row.names = FALSE)

cat("\nMean power by projection (alternatives only):\n")
agg <- aggregate(power ~ method + projection, data = alt, FUN = mean)
print(format(agg[order(agg$method, -agg$power), ], digits = 3), row.names = FALSE)
