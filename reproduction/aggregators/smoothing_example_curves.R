#!/usr/bin/env Rscript
# Curve-level example export backing the section-07 smoothing figures in
# the audit dossier visualisation: three input curves and, for
# each, the two K = 4 polygonal maps the audit contrasts.
#
#   raw       - the 168-point input curve (what the shipped application
#               consumed directly, having lost its reconstruction step).
#   fs_interp - Faber-Schauder truncation: polygonal interpolation at the
#               dyadic vertices k/16 (what the thesis text describes).
#   fs_ls     - least squares onto the same 17-dim polygonal space via
#               Data2fd (what the thesis code did; the corrected app path).
#
# Curves: two ISO weeks of Lithuanian balancing prices selected from the FULL
# fixed-pipeline weekly matrix. The analysis scripts' outlier trim (DROP_WEEK)
# is deliberately NOT applied here, so the removed outlier week stays
# exportable; do not "fix" that. The third curve is the dissertation's
# exemplary spike function
#   f_alpha(x) = (1 - |x - 1/2|^alpha)^(1/alpha),  alpha = 0.45
# (ch1_introduction/subsec_faber_schauder.tex, Fig. fig:three_side_by_side),
# sampled on the same 168-point grid so all three curves run through the
# identical pipeline.
#
# Outputs: data/smoothing_example_curves.csv (committed, full precision)
# and a paste-ready "smooth_example" JSON fragment for the page's DATA
# literal, printed and also written to results/ for copying.
#
# Usage:  Rscript reproduction/smoothing_example_curves.R
# Env:    LITGRID_RAW_DIR overrides the raw-data lookup (default:
#         script until it finds a .raw sibling).

# ---- Inputs (edit these) --------------------------------------------------
WEEKS       <- c(201827L, 202149L)  # ISO year-week labels to export
ROLES       <- c("detected segment start", "removed outlier week")
SPIKE_ALPHA <- 0.45                 # dissertation spike example f_alpha
K           <- 4L                   # dyadic level; space dim = 2^K + 1 = 17
N_INPUT     <- 168L                 # hourly samples per week (fixed by data)
ROUND_WEEK  <- 2L                   # JS decimals for price curves (EUR/MWh)
ROUND_SPIKE <- 4L                   # JS decimals for the spike curve ([0,1])
GAP_SD_BAND <- c(0.2, 2)            # warn if a week's gap/sd leaves this band
OUT_CSV     <- "smoothing_example_curves.csv"        # under data/
OUT_JS      <- "smoothing_example_fragment.js"       # under results/
# ---------------------------------------------------------------------------

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- sub(key, "", args[grepl(key, args)])
  if (length(hit)) return(normalizePath(hit[1]))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))

source(file.path(root, "R", "legacy", "reconstruction.R"))
source(file.path(root, "applications", "litgrid", "data.R"))

litgrid_dir <- Sys.getenv("LITGRID_RAW_DIR", unset = "")
if (!nzchar(litgrid_dir)) {
  litgrid_dir <- normalizePath(file.path(root, "..", "data-raw", "litgrid-raw"),
                               mustWork = FALSE)
}
if (!dir.exists(litgrid_dir)) {
  stop("Litgrid raw CSVs not found at ", litgrid_dir,
       "; run Rscript applications/litgrid/download_litgrid.R or set LITGRID_RAW_DIR.")
}
litgrid_dir <- litgrid_dir

# ---- Assemble the three input curves (columns) ------------------------------
m <- suppressWarnings(load_litgrid_balancing_prices(litgrid_dir))
missing <- setdiff(as.character(WEEKS), rownames(m))
if (length(missing)) {
  stop("Weeks not in the weekly matrix: ", paste(missing, collapse = ", "),
       " (available ", rownames(m)[1], "..", rownames(m)[nrow(m)], ")")
}
stopifnot(ncol(m) == N_INPUT)

t_in  <- seq(0, 1, length.out = N_INPUT)
spike <- (1 - abs(t_in - 0.5)^SPIKE_ALPHA)^(1 / SPIKE_ALPHA)

X <- cbind(t(m[as.character(WEEKS), , drop = FALSE]), f_alpha = spike)

labels <- c(as.character(WEEKS), "f_alpha")
isos   <- c(sprintf("%d-W%02d", WEEKS %/% 100L, WEEKS %% 100L),
            sprintf("f_alpha, alpha = %.2f", SPIKE_ALPHA))
roles  <- c(ROLES, "dissertation spike example")
kinds  <- c(rep("week", length(WEEKS)), "synthetic")
units  <- c(rep("EUR/MWh", length(WEEKS)), "f(t)")
digs   <- c(rep(ROUND_WEEK, length(WEEKS)), ROUND_SPIKE)

# ---- The two maps ------------------------------------------------------------
Vi <- reconstruct(X, "fs_interp", K, grid = "vertices")$X
Vl <- reconstruct(X, "fs_ls",     K, grid = "vertices")$X
t_v <- seq(0, 1, by = 2^(-K))

stopifnot(nrow(Vi) == 2^K + 1L, nrow(Vl) == 2^K + 1L,
          !anyNA(X), !anyNA(Vi), !anyNA(Vl),
          max(abs(Vi - Vl)) > 0)
# fs_interp at t = 0, 1 must coincide with the raw endpoints exactly.
stopifnot(isTRUE(all.equal(as.numeric(Vi[1L, ]), as.numeric(X[1L, ]))),
          isTRUE(all.equal(as.numeric(Vi[2^K + 1L, ]), as.numeric(X[N_INPUT, ]))))

# ---- Diagnostics -------------------------------------------------------------
sd_raw <- apply(X, 2L, sd)
gap    <- apply(abs(Vi - Vl), 2L, max)
gap_sd <- gap / sd_raw
cat(sprintf("%-10s %-28s %12s %16s %8s\n",
            "curve", "role", "sd(raw)", "max|interp-ls|", "gap/sd"))
for (j in seq_along(labels)) {
  cat(sprintf("%-10s %-28s %12.4f %16.4f %8.3f\n",
              labels[j], roles[j], sd_raw[j], gap[j], gap_sd[j]))
  if (kinds[j] == "week" &&
      (gap_sd[j] < GAP_SD_BAND[1] || gap_sd[j] > GAP_SD_BAND[2])) {
    warning(sprintf("gap/sd for week %s is %.3f, outside [%.1f, %.1f]; the note's observed band is 0.39-1.32",
                    labels[j], gap_sd[j], GAP_SD_BAND[1], GAP_SD_BAND[2]))
  }
}

# ---- Tidy CSV (full precision) -----------------------------------------------
curve_rows <- function(j) {
  base <- data.frame(label = labels[j], iso_week = isos[j], role = roles[j],
                     kind = kinds[j])
  rbind(
    cbind(base, series = "raw",       K = NA_integer_, grid = "input",
          idx = seq_len(N_INPUT), t = t_in, value = X[, j]),
    cbind(base, series = "fs_interp", K = K, grid = "vertices",
          idx = seq_along(t_v), t = t_v, value = Vi[, j]),
    cbind(base, series = "fs_ls",     K = K, grid = "vertices",
          idx = seq_along(t_v), t = t_v, value = Vl[, j]))
}
out <- do.call(rbind, lapply(seq_along(labels), curve_rows))
out_path <- file.path(root, "data", OUT_CSV)
write.csv(out, out_path, row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %s (%d data rows)\n", out_path, nrow(out)))

# ---- Paste-ready JS fragment for the page's DATA literal ---------------------
vec_json <- function(x, d) {
  paste0("[", paste(sprintf("%.*f", d, x), collapse = ","), "]")
}
curve_json <- function(j) {
  sprintf(paste0('{"label":"%s","iso":"%s","role":"%s","kind":"%s",',
                 '"unit":"%s","sd_raw":%s,"gap_sd":%.3f,',
                 '"raw":%s,"fs_interp":%s,"fs_ls":%s}'),
          labels[j], isos[j], roles[j], kinds[j], units[j],
          sprintf("%.*f", digs[j], sd_raw[j]), gap_sd[j],
          vec_json(X[, j],  digs[j]),
          vec_json(Vi[, j], digs[j]),
          vec_json(Vl[, j], digs[j]))
}
fragment <- sprintf(', "smooth_example": {"K": %d, "curves": [%s]}',
                    K, paste(vapply(seq_along(labels), curve_json, ""),
                             collapse = ", "))

js_path <- file.path(root, "results", OUT_JS)
dir.create(dirname(js_path), showWarnings = FALSE, recursive = TRUE)
writeLines(fragment, js_path)
cat(sprintf("Wrote %s (%d chars)\n", js_path, nchar(fragment)))
cat("\n---- paste into DATA, before the final `};` on the const DATA line ----\n")
cat(fragment)
cat("\n---- end paste ----------------------------------------------------------\n")
