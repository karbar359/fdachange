#!/usr/bin/env Rscript
# Turn data/smoothing_sensitivity.csv into the LaTeX fragment included by
# notes/smoothing_reconstruction/main.tex. Follows the same pattern as the
# other aggregate_*_tex.R scripts: read the committed tidy CSV, emit a table,
# never hand-transcribe numbers into the document.
#
# Usage:  Rscript reproduction/aggregate_smoothing_sensitivity_tex.R

# ---- Inputs (edit these) --------------------------------------------------
IN_CSV  <- "data/smoothing_sensitivity.csv"
OUT_TEX <- "output/smoothing_results_table.tex"
# The deterministic projection carries the cleanest signal; the random ones are
# summarised in prose because their redraw variability dominates.
MAIN_PROJECTION <- "l2"
GAMMA_MAIN      <- 0.05
# ---------------------------------------------------------------------------

root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) dirname(dirname(normalizePath(hit[1]))) else getwd()
})

d <- read.csv(file.path(root, IN_CSV), stringsAsFactors = FALSE)

# Order the reconstructions from coarsest to finest, raw last.
recon_order <- function(x) {
  meth <- sub("_K.*$", "", x)
  K    <- suppressWarnings(as.integer(sub("^.*_K([0-9]+)_.*$", "\\1", x)))
  grid <- sub("^.*_", "", x)
  rank_meth <- match(meth, c("fs_interp", "fs_ls", "bspline_pen", "none"))
  order(ifelse(is.na(K), 99L, K), rank_meth, grid, na.last = TRUE)
}

app <- subset(d, arm == "application" & projection == MAIN_PROJECTION)
app <- app[recon_order(app$recon), ]
wide <- reshape(app[, c("recon", "gamma", "stat", "reject_rate", "lo", "hi")],
                timevar = "gamma", idvar = "recon", direction = "wide")

esc <- function(s) gsub("_", "\\\\_", s)

lines <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Litgrid application, Wilcoxon scan with the $\\|\\cdot\\|_{L_2}$ projection on the deduplicated series, across the reconstruction sweep. $T$ is the statistic, $[\\hat\\ell_-,\\hat\\ell_+]$ the detected window as fractions of $n=344$. Critical values are 1.867, 2.538 and 3.851. Reconstructions are ordered from coarsest to finest, with the raw sampling grid last.}",
  "\\label{tab:smoothing-sensitivity}",
  "\\footnotesize",
  "\\begin{tabular}{lrrrcc}",
  "\\toprule",
  "\\textbf{Reconstruction} & \\boldmath$T_{0.05}$ & \\boldmath$T_{0.25}$ & \\boldmath$T_{0.45}$ & \\boldmath$\\hat\\ell_-$ & \\boldmath$\\hat\\ell_+$ \\\\",
  "\\midrule")

for (i in seq_len(nrow(wide))) {
  r <- wide[i, ]
  lines <- c(lines, sprintf(
    "%s & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
    esc(r$recon), r$`stat.0.05`, r$`stat.0.25`, r$`stat.0.45`,
    r$`lo.0.05`, r$`hi.0.05`))
}

lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}", "")

# ---- simulation arm -------------------------------------------------------
sim <- subset(d, arm %in% c("simulation_size", "simulation_power") &
                 gamma == GAMMA_MAIN)
if (nrow(sim)) {
  sw <- reshape(sim[, c("recon", "a", "reject_rate")],
                timevar = "a", idvar = "recon", direction = "wide")
  sw <- sw[recon_order(sw$recon), ]
  acols <- grep("^reject_rate\\.", names(sw), value = TRUE)
  avals <- sub("^reject_rate\\.", "", acols)
  lines <- c(lines,
    "\\begin{table}[!ht]", "\\centering",
    sprintf("\\caption{Scenario F1, Wilcoxon scan with the $\\|\\cdot\\|_{L_2}$ projection at $\\gamma=%.2f$, $n=%d$, %d replicates. Column $a=0$ is empirical size, the others are power. The reconstruction choice is immaterial here.}",
            GAMMA_MAIN, sim$n[1], sim$nrep[1]),
    "\\label{tab:smoothing-sensitivity-sim}", "\\footnotesize",
    sprintf("\\begin{tabular}{l%s}", paste(rep("r", length(acols)), collapse = "")),
    "\\toprule",
    sprintf("\\textbf{Reconstruction} & %s \\\\",
            paste(sprintf("\\boldmath$a=%s$", avals), collapse = " & ")),
    "\\midrule")
  for (i in seq_len(nrow(sw))) {
    vals <- vapply(acols, function(cc) sprintf("%.3f", sw[i, cc]), character(1))
    lines <- c(lines, sprintf("%s & %s \\\\", esc(sw$recon[i]),
                              paste(vals, collapse = " & ")))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

out_path <- normalizePath(file.path(root, OUT_TEX), mustWork = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writeLines(lines, out_path)
cat(sprintf("wrote %s (%d application rows, %d simulation rows)\n",
            out_path, nrow(wide), if (exists("sw")) nrow(sw) else 0L))
