#!/usr/bin/env Rscript
# Generate the dissertation table fragment
# `output/tables/_stability_summary.tex`
# from `results/stability_summary.csv` (produced by code/stability.R).
#
# Each row of the input CSV is one (scenario, a, k_sample, gamma) cell with
# the within-X rejection rate R_hat over N_PROJ Ito draws on a fixed X. The
# table aggregates across the K independently drawn fixed samples per
# (scenario, a, gamma): we report mean(R_hat) and sd(R_hat) across the K
# rows. The reported sd mixes the within-X binomial component
# (sqrt(p(1-p)/N_PROJ)) and the between-X spread; on most cells the latter
# dominates.
#
# Output table style matches `_f4_rkhs_quantiles.tex`: booktabs, footnotesize,
# tight tabcolsep, \multirow grouping. Label `tab:stability:summary`.
#
# Usage:
#   Rscript reproduction/aggregate_stability.R

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
out_dir <- normalizePath(
  file.path(root, "output", "tables"))

src <- file.path(res_dir, "stability_summary.csv")
if (!file.exists(src)) stop("Input CSV not found: ", src)
df <- read.csv(src, stringsAsFactors = FALSE)

SCENARIOS <- c("F1", "F2", "F3", "F4")
df <- df[df$scenario %in% SCENARIOS, , drop = FALSE]
if (!nrow(df)) stop("No rows for scenarios ", paste(SCENARIOS, collapse = ","))

n_val      <- max(df$n)
n_proj_val <- max(df$n_proj)
K_val      <- length(unique(df$k_sample))
gamma_set  <- sort(unique(df$gamma))

fmt_a   <- function(a) {
  if (abs(a - round(a)) < 1e-9) sprintf("%d", as.integer(round(a)))
  else                          sprintf("%g", a)
}
fmt_g   <- function(g) if (g == 0) "0" else sprintf("%g", g)
fmt_R   <- function(x) if (!is.finite(x)) "--" else sprintf("%.3f", x)

emit_rows <- function() {
  body <- character(0)
  for (si in seq_along(SCENARIOS)) {
    sc <- SCENARIOS[si]
    sub_sc <- df[df$scenario == sc, , drop = FALSE]
    if (!nrow(sub_sc)) next
    a_vals <- sort(unique(sub_sc$a))
    n_a    <- length(a_vals)
    n_g    <- length(gamma_set)
    rows_per_sc <- n_a * n_g
    first_in_sc <- TRUE
    for (a in a_vals) {
      sub_a <- sub_sc[abs(sub_sc$a - a) < 1e-9, , drop = FALSE]
      first_in_a <- TRUE
      for (g in gamma_set) {
        cell <- sub_a[abs(sub_a$gamma - g) < 1e-9, , drop = FALSE]
        if (nrow(cell)) {
          R_mean <- mean(cell$R_hat)
          R_sd   <- if (nrow(cell) > 1L) sd(cell$R_hat) else 0
        } else {
          R_mean <- NA_real_; R_sd <- NA_real_
        }
        left  <- if (first_in_sc) sprintf("\\multirow{%d}{*}{%s}", rows_per_sc, sc) else ""
        amid  <- if (first_in_a)  sprintf("\\multirow{%d}{*}{%s}", n_g, fmt_a(a))  else ""
        body <- c(body, sprintf(
          "%s & %s & %s & %s & %s \\\\",
          left, amid, fmt_g(g), fmt_R(R_mean), fmt_R(R_sd)))
        first_in_sc <- FALSE
        first_in_a  <- FALSE
      }
    }
    if (si < length(SCENARIOS)) body <- c(body, "\\midrule")
  }
  body
}

cap <- sprintf(paste0(
  "\\caption{Stability of the antisymmetric-kernel Wilcoxon scan under ",
  "repeated draws of the It\\^o projector. For each scenario and value of ",
  "$a$, $K=%d$ independently drawn samples of size $n=%d$ are fixed and ",
  "$N=%d$ It\\^o trajectories are drawn against each; the table reports ",
  "the mean and standard deviation of $\\widehat R_N=N^{-1}\\sum_i\\mathbf 1",
  "\\{T_n(\\omega_i)>c_\\alpha\\}$ across the $K$ fixed samples at ",
  "$\\alpha=0.05$.}\\label{tab:stability:summary}"),
  K_val, n_val, n_proj_val)

header <- c(
  "Scenario & $a$ & $\\gamma$ & mean $\\widehat R_N$ & sd $\\widehat R_N$ \\\\")

tex <- c(
  "\\begin{table}[ht]", "\\centering", "\\footnotesize",
  "\\setlength{\\tabcolsep}{4pt}", cap,
  "\\begin{tabular}{l r r r r}", "\\toprule",
  header, "\\midrule",
  emit_rows(),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")

out_path <- file.path(out_dir, "_stability_summary.tex")
writeLines(tex, out_path)
cat(sprintf("wrote %s  (K=%d, n=%d, N_proj=%d)\n",
            out_path, K_val, n_val, n_proj_val))
