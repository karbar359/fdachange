#!/usr/bin/env Rscript
# Generate dissertation table fragments _pareto_F<i>_<method>_<short|long>.tex
# from per-cell driver CSVs in results/.
#
# Reads (lp|wilcoxon|rkhs)_F[1-3]*.csv from results/, filters to n = N
# (default 400), and writes one .tex fragment per (scenario, method,
# segment_length).
#
# Usage:
#   Rscript reproduction/aggregate_pareto_curves_tex.R

# ---- Inputs (edit these to retune) -----------------------------------------
# Pre-set any of these in .GlobalEnv before sourcing to retarget the run; the
# corrected-scan tables are produced by setting SCAN_FAMILY <- "exact" and
# OUT_DIR somewhere else if needed, since the printed chapter is
# frozen and must keep matching the print.
if (!exists("N_TARGET"))   N_TARGET   <- 400L   # only rows with this n end up in the tables
if (!exists("SCAN_FAMILY")) SCAN_FAMILY <- "legacy"  # legacy | exact | corrected
if (!exists("OUT_DIR"))    OUT_DIR    <- NULL   # NULL -> output/tables
GAMMAS   <- c(0, 0.25, 0.45)
WILCOXON_M_LEVELS <- c(1, 2, 4)
A_GRID_BY_SCENARIO <- list(
  F1 = c(1.5, 2, 2.5, 3, 5, 10),
  F2 = c(0.5, 1, 1.5, 2, 5),
  F3 = c(0.5, 1, 1.5, 2, 2.5, 5)
)
IDENTITY <- list(
  F1 = "Normal vs centred Pareto",
  F2 = "Pareto vs Pareto",
  F3 = "Normal vs heavy-tailed (piecewise)"
)
# ---------------------------------------------------------------------------

PROJ_LABEL <- list(
  ito       = "It\\^o",
  l2        = "$\\|\\cdot\\|_{L_2}$",
  c_max_var = "$C$ max-var"
)
PROJ_ORDER <- c("ito", "l2", "c_max_var")
KERNEL_LABEL <- list(
  gaussian = "Gaussian",
  min      = "min",
  `2_max`  = "$2-\\max$"
)
RKHS_PAIRS <- list(
  list(projection = "ito",       kernel = "gaussian"),
  list(projection = "l2",        kernel = "min"),
  list(projection = "c_max_var", kernel = "gaussian")
)

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
out_dir <- if (is.null(OUT_DIR)) {
  normalizePath(file.path(root, "output", "tables"))
} else {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  normalizePath(OUT_DIR)
}

files <- list.files(res_dir,
                    pattern = "^(lp|wilcoxon|rkhs)_F[1-3].*\\.csv$",
                    full.names = TRUE)
files <- files[!grepl("preswap_backup", files)]
# `_exact` CSVs are the corrected-scan re-run (run_exact_scan_n400.sh). The two
# families share every join key, so rbind order would decide which value wins;
# they must never be read together. SCAN_FAMILY picks exactly one. Compare them
# with compare_exact_vs_legacy.R instead of merging.
# `legacybase` CSVs exist only to give compare_exact_vs_legacy.R a faithful
# baseline for the two cells that results/ has since overwritten; never input.
files <- files[!grepl("legacybase", basename(files))]

bn <- basename(files)
files <- switch(SCAN_FAMILY,
  # Only Wilcoxon has an exact re-run; Lp and RKHS do not use the scan and
  # carry over unchanged.
  exact     = files[grepl("_exact", bn) |
                    (!grepl("^wilcoxon_", bn) & !grepl("_corrected", bn))],
  # The corrected family exists for every method.
  corrected = files[grepl("_corrected", bn)],
  legacy    = files[!grepl("_exact|_corrected", bn)],
  stop("SCAN_FAMILY must be legacy | exact | corrected, got ", SCAN_FAMILY))
if (!length(files)) stop("No F1/F2/F3 result CSVs in ", res_dir)
cat(sprintf("Reading %d CSV(s) from %s\n", length(files), res_dir))
long <- do.call(rbind, lapply(files, function(f) {
  read.csv(f, stringsAsFactors = FALSE)
}))
long <- long[long$n == N_TARGET, , drop = FALSE]
if (!nrow(long)) stop("No rows with n=", N_TARGET, " in result CSVs.")

format_pow <- function(p) if (is.na(p)) "--" else sprintf("%.2f", p)
gamma_str  <- function(g) if (g == 0) "0" else sprintf("%g", g)

pick_powers <- function(rows_df, a_grid) {
  vapply(a_grid, function(av) {
    r <- rows_df[abs(rows_df$a - av) < 1e-9, ]
    if (nrow(r)) r$power[1] else NA_real_
  }, numeric(1))
}

emit_wilcoxon <- function(scenario, segment_label, identity_name, rows_df,
                          a_grid, n_val, n_rep, alpha) {
  ell  <- if (segment_label == "short") "\\lfloor\\log n\\rfloor+1" else "n/5"
  cap  <- sprintf(
    "\\caption{Wilcoxon test power on \\emph{%s: %s}, segment length $\\ell^*=%s$. $N_{\\rm rep}=%d$, $n=%d$, $\\alpha=%g$.}\\label{tab:pareto-%s-wilcoxon-%s}",
    scenario, identity_name, ell, n_rep, n_val, alpha,
    tolower(scenario), segment_label)
  align <- paste(c(rep("l", 3), rep("r", length(a_grid))), collapse = " ")
  header <- sprintf("Projection & $m$ & $\\gamma$ & %s \\\\",
                    paste(sprintf("$a=%s$", a_grid), collapse = " & "))
  body <- character(0)
  for (pi_ in seq_along(PROJ_ORDER)) {
    proj <- PROJ_ORDER[pi_]
    first_proj <- TRUE
    for (m in WILCOXON_M_LEVELS) {
      first_m <- TRUE
      for (gamma in GAMMAS) {
        sub <- rows_df[rows_df$projection == proj & rows_df$m == m &
                       abs(rows_df$gamma - gamma) < 1e-9, , drop = FALSE]
        powers <- pick_powers(sub, a_grid)
        left <- if (first_proj) sprintf("\\multirow{9}{*}{%s}", PROJ_LABEL[[proj]]) else ""
        mid  <- if (first_m) sprintf("\\multirow{3}{*}{%d}", m) else ""
        body <- c(body, sprintf("%s & %s & %s & %s \\\\",
                                left, mid, gamma_str(gamma),
                                paste(vapply(powers, format_pow, character(1)),
                                      collapse = " & ")))
        first_proj <- FALSE
        first_m    <- FALSE
      }
    }
    if (pi_ < length(PROJ_ORDER)) body <- c(body, "\\midrule")
  }
  c("\\begin{table}[!ht]", "\\centering", "\\footnotesize",
    "\\setlength{\\tabcolsep}{3.5pt}", cap,
    sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule",
    body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

emit_lp <- function(scenario, segment_label, identity_name, rows_df,
                    a_grid, n_val, n_rep, alpha) {
  ell  <- if (segment_label == "short") "\\lfloor\\log n\\rfloor+1" else "n/5"
  cap  <- sprintf(
    "\\caption{$L_p$ test power on \\emph{%s: %s}, segment length $\\ell^*=%s$. $N_{\\rm rep}=%d$, $n=%d$, $\\alpha=%g$.}\\label{tab:pareto-%s-lp-%s}",
    scenario, identity_name, ell, n_rep, n_val, alpha,
    tolower(scenario), segment_label)
  align <- paste(c(rep("l", 3), rep("r", length(a_grid))), collapse = " ")
  header <- sprintf("Projection & $p$ & $\\gamma$ & %s \\\\",
                    paste(sprintf("$a=%s$", a_grid), collapse = " & "))
  body <- character(0)
  for (pi_ in seq_along(PROJ_ORDER)) {
    proj <- PROJ_ORDER[pi_]
    first_proj <- TRUE
    for (gamma in GAMMAS) {
      sub <- rows_df[rows_df$projection == proj &
                     abs(rows_df$gamma - gamma) < 1e-9, , drop = FALSE]
      powers <- pick_powers(sub, a_grid)
      left <- if (first_proj) sprintf("\\multirow{3}{*}{%s}", PROJ_LABEL[[proj]]) else ""
      mid  <- if (first_proj) "\\multirow{3}{*}{2}" else ""
      body <- c(body, sprintf("%s & %s & %s & %s \\\\",
                              left, mid, gamma_str(gamma),
                              paste(vapply(powers, format_pow, character(1)),
                                    collapse = " & ")))
      first_proj <- FALSE
    }
    if (pi_ < length(PROJ_ORDER)) body <- c(body, "\\midrule")
  }
  c("\\begin{table}[!ht]", "\\centering", "\\footnotesize",
    "\\setlength{\\tabcolsep}{3.5pt}", cap,
    sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule",
    body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

emit_rkhs <- function(scenario, segment_label, identity_name, rows_df,
                      a_grid, n_val, n_rep, alpha) {
  ell  <- if (segment_label == "short") "\\lfloor\\log n\\rfloor+1" else "n/5"
  cap  <- sprintf(
    "\\caption{RKHS test power on \\emph{%s: %s}, segment length $\\ell^*=%s$. $N_{\\rm rep}=%d$, $n=%d$, $\\alpha=%g$.}\\label{tab:pareto-%s-rkhs-%s}",
    scenario, identity_name, ell, n_rep, n_val, alpha,
    tolower(scenario), segment_label)
  align <- paste(c(rep("l", 3), rep("r", length(a_grid))), collapse = " ")
  header <- sprintf("Projection & Kernel & $\\gamma$ & %s \\\\",
                    paste(sprintf("$a=%s$", a_grid), collapse = " & "))
  body <- character(0)
  for (pi_ in seq_along(RKHS_PAIRS)) {
    pair <- RKHS_PAIRS[[pi_]]
    first_proj <- TRUE
    for (gamma in GAMMAS) {
      sub <- rows_df[rows_df$projection == pair$projection &
                     abs(rows_df$gamma - gamma) < 1e-9, , drop = FALSE]
      powers <- pick_powers(sub, a_grid)
      left <- if (first_proj) sprintf("\\multirow{3}{*}{%s}", PROJ_LABEL[[pair$projection]]) else ""
      mid  <- if (first_proj) sprintf("\\multirow{3}{*}{%s}", KERNEL_LABEL[[pair$kernel]]) else ""
      body <- c(body, sprintf("%s & %s & %s & %s \\\\",
                              left, mid, gamma_str(gamma),
                              paste(vapply(powers, format_pow, character(1)),
                                    collapse = " & ")))
      first_proj <- FALSE
    }
    if (pi_ < length(RKHS_PAIRS)) body <- c(body, "\\midrule")
  }
  c("\\begin{table}[!ht]", "\\centering", "\\footnotesize",
    "\\setlength{\\tabcolsep}{3.5pt}", cap,
    sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule",
    body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

EMITTERS <- list(wilcoxon = emit_wilcoxon, lp = emit_lp, rkhs = emit_rkhs)

for (scenario in names(A_GRID_BY_SCENARIO)) {
  scen <- long[long$scenario == scenario, , drop = FALSE]
  if (!nrow(scen)) {
    cat(sprintf("[%s] no rows at n=%d, skipping\n", scenario, N_TARGET))
    next
  }
  a_grid <- A_GRID_BY_SCENARIO[[scenario]]
  n_val  <- max(scen$n)
  n_rep  <- max(scen$N_rep)
  alpha  <- max(scen$alpha)
  ec_lens <- sort(unique(scen$ec_length))
  if (length(ec_lens) < 2) {
    cat(sprintf("[%s] only %d ec_length values found; expected 2\n",
                scenario, length(ec_lens)))
  }
  short_ec <- ec_lens[1]
  long_ec  <- ec_lens[length(ec_lens)]
  for (method in names(EMITTERS)) {
    method_rows <- scen[scen$method == method, , drop = FALSE]
    if (!nrow(method_rows)) {
      cat(sprintf("[%s/%s] no rows, skipping\n", scenario, method))
      next
    }
    short_rows <- method_rows[abs(method_rows$ec_length - short_ec) < 1e-6, , drop = FALSE]
    long_rows  <- method_rows[abs(method_rows$ec_length - long_ec)  < 1e-6, , drop = FALSE]
    out_short  <- file.path(out_dir, sprintf("_pareto_%s_%s_short.tex", scenario, method))
    out_long   <- file.path(out_dir, sprintf("_pareto_%s_%s_long.tex",  scenario, method))
    writeLines(EMITTERS[[method]](scenario, "short", IDENTITY[[scenario]],
                                  short_rows, a_grid, n_val, n_rep, alpha), out_short)
    writeLines(EMITTERS[[method]](scenario, "long",  IDENTITY[[scenario]],
                                  long_rows, a_grid, n_val, n_rep, alpha), out_long)
    cat(sprintf("[%s/%s] wrote %s and %s\n", scenario, method,
                basename(out_short), basename(out_long)))
  }
}
cat("Done.\n")
