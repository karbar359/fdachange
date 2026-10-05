#!/usr/bin/env Rscript
# Combine the Lp exponent-study CSVs written by sweeps/run_lp_p_study.sh
# into the committed data/lp_p_power.csv, and print the comparison the
# study exists for: does the exponent p matter for the corrected Lp scan
# on the scalar scenarios?
#
# Usage:  Rscript reproduction/aggregators/aggregate_lp_p.R

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) return(normalizePath(hit[1]))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))

files <- list.files(file.path(root, "results"),
                    pattern = "^lp_S[1-5]_n400_p(2|4|10)_corrected\\.csv$",
                    full.names = TRUE)
if (length(files) < 15L) {
  stop("Expected 15 result CSVs in results/, found ", length(files),
       "; run sweeps/run_lp_p_study.sh first.")
}

out <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
out$p <- as.integer(sub("^p([0-9]+)_corrected$", "\\1", out$tag))
stopifnot(!anyNA(out$p), all(out$p %in% c(2L, 4L, 10L)))
out <- out[order(out$scenario, out$p, out$ec_length, out$gamma, out$a), ]

out_path <- file.path(root, "data", "lp_p_power.csv")
write.csv(out, out_path, row.names = FALSE)
cat(sprintf("wrote %s (%d rows from %d files)\n\n", out_path, nrow(out),
            length(files)))

# ---- the comparison of record ----------------------------------------------
cat("Mean power over a > 0 by (scenario, window, gamma) across p:\n")
alt <- out[out$a > 0, ]
agg <- aggregate(power ~ scenario + ec_length + gamma + p, alt, mean)
wide <- reshape(agg, idvar = c("scenario", "ec_length", "gamma"),
                timevar = "p", direction = "wide")
names(wide) <- sub("^power\\.", "p", names(wide))
wide <- wide[order(wide$scenario, wide$ec_length, wide$gamma), ]
print(format(wide, digits = 2), row.names = FALSE)

cat("\nEmpirical size (a = 0), worst cell per p:\n")
sz <- aggregate(power ~ p, out[out$a == 0, ], max)
print(format(sz, digits = 2), row.names = FALSE)
