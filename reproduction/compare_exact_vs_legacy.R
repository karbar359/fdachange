#!/usr/bin/env Rscript
# Pair a re-run family of simulation results (exact or corrected) with the
# legacy CSVs that fed the printed tables, and report what moved.
#
# FAMILY = "exact":     the Wilcoxon-only re-run that fixes the summed-area
#                       off-by-one. Same statistic family otherwise,
#                       per-cell reseeding, so deltas are exact.
# FAMILY = "corrected": the full audit family across all four methods
#                       (weight direction + index, Lp centring, RKHS printed
#                       formula + regenerated null, Birbilas printed
#                       denominator + regenerated finite null, Kuiper gamma=0
#                       null, exact-length EC windows). The EC window change
#                       alters the draw counts, so the RNG streams shift:
#                       corrected-vs-legacy deltas here include MC noise of
#                       order 1/sqrt(N_rep) on top of the statistic change.
#
# Pairing is explicit, never by glob: for each variant CSV we look for a
# legacy CSV with the same (method, scenario, n, m, projection, tag-stem,
# N_rep) whose (gamma, ec_length, a) grid matches exactly. Unpaired files
# are reported rather than silently dropped. N_rep is part of the key on
# purpose: `results/` is working state and cells can be overwritten by
# higher-N_rep reruns.
#
# The committed reference outputs live in data/{exact,corrected}_vs_legacy.csv;
# a full re-run of the sweeps followed by this script must reproduce them.
#
# Usage:  Rscript reproduction/compare_exact_vs_legacy.R              (exact)
#         Rscript -e "FAMILY <- 'corrected'; source('reproduction/compare_exact_vs_legacy.R')"

# ---- Inputs (edit these) --------------------------------------------------
if (!exists("FAMILY"))   FAMILY   <- "exact"    # exact | corrected
RESULTS_DIR <- "results"
N_TARGET    <- 400L
# Power is printed to 2 decimals, so a cell "changes in print" if it rounds
# differently.
PRINT_DIGITS <- 2L
# ---------------------------------------------------------------------------

stopifnot(FAMILY %in% c("exact", "corrected"))
OUT_CSV <- sprintf("data/%s_vs_legacy.csv", FAMILY)
ST_CSV  <- sprintf("data/%s_vs_legacy_stability.csv", FAMILY)

root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) dirname(normalizePath(hit[1])) else getwd()
})
res_dir <- file.path(root, RESULTS_DIR)

key_cols  <- c("gamma", "ec_length", "a")
meta_cols <- c("method", "scenario", "n", "m", "projection")

read_cell <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$file <- basename(path)
  d
}

prefix <- if (FAMILY == "exact") "^wilcoxon_" else "^(wilcoxon|lp|rkhs|birbilas)_"
all_files <- list.files(res_dir, pattern = paste0(prefix, ".*\\.csv$"),
                        full.names = TRUE)
all_files <- all_files[!grepl("preswap_backup", all_files)]
tagpat <- paste0("_", FAMILY)
variant_files <- grep(tagpat, all_files, value = TRUE, fixed = TRUE)
legacy_files  <- all_files[!grepl("_exact|_corrected", basename(all_files))]

if (!length(variant_files)) stop("No *", tagpat, "* CSVs in ", res_dir)

legacy <- lapply(legacy_files, read_cell)
names(legacy) <- basename(legacy_files)

grid_of <- function(d) paste(sort(apply(d[, key_cols], 1, paste, collapse = "|")),
                             collapse = ";")
# The tag with the family marker stripped: pairs e.g. "p2_corrected" <-> "p2",
# "min_runif_l2_corrected" <-> "min_runif_l2", "h_0_beta0_exact" <-> "h_0_beta0".
stem_of <- function(tag) gsub("_(exact|corrected)", "", tag)

rows <- list()
unpaired <- character(0)

for (f in sort(variant_files)) {
  e <- read_cell(f)
  e <- e[e$n == N_TARGET, , drop = FALSE]
  if (!nrow(e)) next
  cands <- Filter(function(d) {
    d <- d[d$n == N_TARGET, , drop = FALSE]
    nrow(d) > 0 &&
      d$method[1] == e$method[1] &&
      d$scenario[1] == e$scenario[1] && d$m[1] == e$m[1] &&
      d$projection[1] == e$projection[1] &&
      stem_of(d$tag[1]) == stem_of(e$tag[1]) &&
      identical(sort(unique(d$N_rep)), sort(unique(e$N_rep))) &&
      grid_of(d) == grid_of(e)
  }, legacy)

  if (length(cands) != 1L) {
    unpaired <- c(unpaired, sprintf("%s -> %d candidates", basename(f), length(cands)))
    next
  }
  l <- cands[[1]]
  l <- l[l$n == N_TARGET, , drop = FALSE]

  mrg <- merge(l[, c(meta_cols, key_cols, "N_rep", "power")],
               e[, c(key_cols, "power")],
               by = key_cols, suffixes = c("_legacy", "_new"))
  mrg$delta <- mrg$power_new - mrg$power_legacy
  mrg$legacy_file <- names(cands)[1]
  mrg$new_file    <- basename(f)
  rows[[length(rows) + 1L]] <- mrg
}

out <- do.call(rbind, rows)
if (is.null(out) || !nrow(out)) stop("No pairings found.")

out <- out[order(out$method, out$scenario, out$projection, out$m,
                 out$ec_length, out$gamma, out$a), ]
out$changes_in_print <- round(out$power_legacy, PRINT_DIGITS) !=
                        round(out$power_new,  PRINT_DIGITS)

out_path <- file.path(root, OUT_CSV)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_path, row.names = FALSE)

# ---- report ---------------------------------------------------------------
cat(sprintf("[%s vs legacy] paired %d cells across %d file pairs.\n",
            FAMILY, nrow(out), length(unique(out$new_file))))
if (length(unpaired)) {
  cat("UNPAIRED (inspect before trusting the summary):\n")
  cat(paste0("  ", unpaired, collapse = "\n"), "\n")
}
if (FAMILY == "corrected") {
  cat("NOTE: EC_EXACT shifts the RNG streams, so per-cell deltas carry MC\n")
  cat("      noise of order 1/sqrt(N_rep) in addition to the statistic change.\n")
}

cat(sprintf("\nPower change: mean delta %+.4f, mean |delta| %.4f, max |delta| %.4f\n",
            mean(out$delta), mean(abs(out$delta)), max(abs(out$delta))))
cat(sprintf("Cells that differ at all: %d of %d (%.1f%%)\n",
            sum(out$delta != 0), nrow(out), 100 * mean(out$delta != 0)))
cat(sprintf("Cells that change at %d decimals (i.e. in print): %d (%.1f%%)\n",
            PRINT_DIGITS, sum(out$changes_in_print),
            100 * mean(out$changes_in_print)))
cat(sprintf("Direction among changed cells: %d up, %d down\n",
            sum(out$delta > 0), sum(out$delta < 0)))

cat("\nBy method:\n")
aggm <- aggregate(cbind(delta, abs_delta = abs(delta)) ~ method, data = out, FUN = mean)
aggm$n_cells <- as.vector(table(out$method)[aggm$method])
aggm$changed <- as.vector(tapply(out$changes_in_print, out$method, sum)[aggm$method])
print(format(aggm, digits = 3), row.names = FALSE)

cat("\nBy scenario:\n")
agg <- aggregate(cbind(delta, abs_delta = abs(delta)) ~ scenario, data = out, FUN = mean)
agg$n_cells <- as.vector(table(out$scenario)[agg$scenario])
agg$changed <- as.vector(tapply(out$changes_in_print, out$scenario, sum)[agg$scenario])
print(format(agg, digits = 3), row.names = FALSE)

cat("\nBy gamma:\n")
agg2 <- aggregate(cbind(delta, abs_delta = abs(delta)) ~ gamma, data = out, FUN = mean)
agg2$changed <- as.vector(tapply(out$changes_in_print, out$gamma, sum)[as.character(agg2$gamma)])
print(format(agg2, digits = 3), row.names = FALSE)

big <- out[abs(out$delta) >= 0.10, ]
if (nrow(big)) {
  cat(sprintf("\nCells moving by >= 0.10 (%d):\n", nrow(big)))
  print(format(big[order(-abs(big$delta)),
                   c("method","scenario","projection","m","gamma","ec_length","a",
                     "power_legacy","power_new","delta")], digits = 3),
        row.names = FALSE)
} else {
  cat("\nNo cell moves by 0.10 or more.\n")
}

cat(sprintf("\nwrote %s\n", out_path))

# ---- stability table ------------------------------------------------------
# The stability study reports the mean and sd of R_hat_N across K fixed
# samples, so it needs its own pairing.
stab_l <- file.path(res_dir, "stability_summary.csv")
stab_e <- file.path(res_dir, sprintf("stability_summary_%s.csv", FAMILY))
if (file.exists(stab_l) && file.exists(stab_e)) {
  agg <- function(f) {
    d <- read.csv(f)
    s <- aggregate(R_hat ~ scenario + a + gamma, data = d,
                   FUN = function(x) c(m = mean(x), s = sd(x)))
    s <- do.call(data.frame, s)
    names(s)[4:5] <- c("mean_R", "sd_R")
    s
  }
  st <- merge(agg(stab_l), agg(stab_e), by = c("scenario", "a", "gamma"),
              suffixes = c("_legacy", "_new"))
  st$delta_mean <- st$mean_R_new - st$mean_R_legacy
  st <- st[order(st$scenario, st$a, st$gamma), ]
  st_path <- file.path(root, ST_CSV)
  write.csv(st, st_path, row.names = FALSE)
  cat(sprintf("\nStability: %d cells, mean |delta| %.4f, max %.4f, %d change at 3 dp\n",
              nrow(st), mean(abs(st$delta_mean)), max(abs(st$delta_mean)),
              sum(round(st$mean_R_legacy, 3) != round(st$mean_R_new, 3))))
  cat(sprintf("wrote %s\n", st_path))
} else {
  cat(sprintf("\nStability comparison skipped (need results/stability_summary_%s.csv)\n",
              FAMILY))
  st <- NULL
}

# ---- LaTeX summary fragment -----------------------------------------------
TEX_OUT <- file.path(root, "output",
                     sprintf("%s_vs_legacy_table.tex", FAMILY))
dir.create(dirname(TEX_OUT), showWarnings = FALSE, recursive = TRUE)
by_m <- aggregate(cbind(delta, abs_delta = abs(delta)) ~ m, data = out, FUN = mean)
by_m$n_cells <- as.vector(table(out$m)[as.character(by_m$m)])
by_m$changed <- as.vector(tapply(out$changes_in_print, out$m, sum)[as.character(by_m$m)])

cap <- if (FAMILY == "exact") {
  sprintf("\\caption{Effect of the corrected scan on the simulation chapter. All %d Wilcoxon power cells were re-run at the printed parameterisation and paired with their legacy counterparts. The bias is $O(1/n)$ in the per-block sample size, so it grows with the split $m$: at $m=4$ each block holds only $n/m=100$ observations. Power at $\\gamma=0.45$ is zero under both scans, so nothing there can move.}", nrow(out))
} else {
  sprintf("\\caption{Effect of the corrected statistics on the simulation chapter: %d cells across all four methods, re-run at the printed parameterisation with the printed formulas, recalibrated nulls, and exact-length change windows, paired with their legacy counterparts. Deltas include Monte-Carlo noise of order $1/\\sqrt{N_{\\rm rep}}$ because the window fix shifts the RNG streams.}", nrow(out))
}
tex <- c(
  "\\begin{table}[!ht]", "\\centering", cap,
  sprintf("\\label{tab:%s-vs-legacy}", FAMILY), "\\footnotesize",
  "\\begin{tabular}{lrrrr}", "\\toprule",
  "\\textbf{Split} & \\textbf{block size} & \\textbf{cells} & \\textbf{mean $\\Delta$ power} & \\textbf{changed in print} \\\\",
  "\\midrule")
for (i in seq_len(nrow(by_m))) {
  tex <- c(tex, sprintf("$m=%d$ & %d & %d & $%+.4f$ & %d (%.1f\\%%) \\\\",
                        by_m$m[i], 400L %/% by_m$m[i], by_m$n_cells[i],
                        by_m$delta[i], by_m$changed[i],
                        100 * by_m$changed[i] / by_m$n_cells[i]))
}
tex <- c(tex, "\\midrule",
  sprintf("\\textbf{all} & & %d & $%+.4f$ & %d (%.1f\\%%) \\\\", nrow(out),
          mean(out$delta), sum(out$changes_in_print),
          100 * mean(out$changes_in_print)),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, TEX_OUT)
cat(sprintf("wrote %s\n", normalizePath(TEX_OUT)))
