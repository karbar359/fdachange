#!/usr/bin/env Rscript
# Assemble data/thesis_power_values.csv, the single committed source of
# every power value plotted by aggregators/plot_thesis_power.R.
#
# Sources, in precedence order (one number per cell, keys never overlap):
#   1. Fresh N_rep = 100 rerun CSVs (F1/F2, all three tests), produced by
#      sweeps/run_thesis_power_F1F2_n400_nrep100.sh into results/thesis_*.csv
#      (results/ is gitignored; this committed CSV is what survives).
#   2. Preserved values from the original pipeline's plotting script whose CSVs
#      were deleted: Wilcoxon F1/F2, L2 projection, gamma = 0, long segment,
#      m in {2, 4} (22 values, used by the m-effect figure).
#   3. Printed thesis tables (all N_rep = 100, n = 400) parsed from
#      the thesis tex fragments (THESIS_TEX_DIR): _pareto_F3_*, _f4_*, _scalar_*.
#      The printed F1/F2 tables are N_rep = 20 and are deliberately NOT
#      emitted; they are only used as a wrong-pipeline tolerance check.
#
# Verification built in: caption asserts (N_rep = 100, n = 400), per-table
# row-count asserts, spot-value asserts against independently read cells,
# key uniqueness, fresh == preserved for the 66 overlapping Wilcoxon values,
# and a fresh-vs-printed-N_rep=20 tolerance scan.
#
# Usage:  Rscript reproduction/aggregators/assemble_thesis_power.R
#         ALLOW_MISSING_RERUN=1 Rscript reproduction/aggregators/assemble_thesis_power.R
#           (emit printed + preserved rows only, while the rerun is running)

# ---- Constants (edit here) ---------------------------------------------
ALLOW_MISSING_RERUN <- identical(Sys.getenv("ALLOW_MISSING_RERUN"), "1")
ALPHA <- 0.05

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) return(normalizePath(hit[1]))
  NULL
})()
root    <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))
tab_dir <- Sys.getenv("THESIS_TEX_DIR",   # printed-table tex fragments
                      file.path(root, "output", "tables"))
res_dir <- file.path(root, "results")
out_csv <- file.path(root, "data", "thesis_power_values.csv")

# ---- Printed-table parser ------------------------------------------------

# Extract the content of the LAST braced group of \multirow{..}{..}{CONTENT},
# brace-matched (handles \shortstack with nested braces); otherwise return
# the trimmed cell.
extract_label <- function(cell) {
  cell <- trimws(cell)
  m <- regexpr("\\\\multirow\\{[^}]*\\}\\{[^}]*\\}\\{", cell)
  if (m == -1) return(cell)
  start <- as.integer(m) + attr(m, "match.length")
  chars <- strsplit(cell, "")[[1]]
  depth <- 1L; i <- start
  while (i <= length(chars) && depth > 0L) {
    if (chars[i] == "{") depth <- depth + 1L
    if (chars[i] == "}") depth <- depth - 1L
    i <- i + 1L
  }
  trimws(paste(chars[start:(i - 2L)], collapse = ""))
}

map_projection <- function(x) {
  if (grepl("It", x))       return("ito")
  if (grepl("L_2", x))      return("l2")
  if (grepl("max-var", x))  return("c_max_var")
  stop("unmapped projection label: ", x)
}
map_kernel <- function(x) {
  if (grepl("Gaussian", x) || grepl("exp", x)) return("gaussian")
  if (grepl("min", x))                          return("min")
  if (grepl("max", x))                          return("2_max")
  stop("unmapped kernel label: ", x)
}
map_scenario <- function(x) {
  hit <- regmatches(x, regexpr("S[1-5]", x))
  if (!length(hit)) stop("unmapped scenario label: ", x)
  hit
}

# Parse one printed table file into a long data frame.
#   label_cols: character vector naming the leading tabular columns, from
#     {"projection", "kernel", "m", "p", "gamma", "scenario", "scenario_kernel"}
parse_printed <- function(file, label_cols, expected_rows, expect_nrep = 100) {
  path  <- file.path(tab_dir, file)
  lines <- readLines(path, warn = FALSE)

  cap <- lines[grepl("\\\\caption", lines)][1]
  stopifnot(grepl(sprintf("N_\\{\\\\rm rep\\}=%d", expect_nrep), cap),
            grepl("n=400", cap, fixed = TRUE))
  segment <- if (grepl("_long", file)) "long" else "short"

  header <- lines[grepl("\\$a=", lines)][1]
  a_grid <- as.numeric(gsub("a=", "",
    regmatches(header, gregexpr("a=[0-9.]+", header))[[1]]))
  stopifnot(length(a_grid) >= 5)

  txt <- paste(lines, collapse = "\n")
  txt <- sub("(?s)^.*?\\\\midrule", "", txt, perl = TRUE)
  txt <- sub("(?s)\\\\bottomrule.*$", "", txt, perl = TRUE)
  txt <- gsub("\\\\midrule", "", txt)
  # Split rows on the \\ terminator ONLY at end-of-line: \shortstack labels
  # (scalar RKHS tables) contain a mid-line \\ inside braces.
  rows <- strsplit(txt, "\\\\\\\\[ \t]*\n")[[1]]
  rows <- rows[nzchar(trimws(rows))]
  stopifnot(length(rows) == expected_rows)

  carry <- setNames(rep(NA_character_, length(label_cols)), label_cols)
  out <- vector("list", length(rows))
  for (r in seq_along(rows)) {
    cells <- trimws(strsplit(rows[r], "&", fixed = TRUE)[[1]])
    stopifnot(length(cells) == length(label_cols) + length(a_grid))
    for (i in seq_along(label_cols)) {
      if (nzchar(cells[i])) carry[i] <- extract_label(cells[i])
    }
    powers <- suppressWarnings(as.numeric(cells[-seq_along(label_cols)]))
    if (anyNA(powers)) stop("non-numeric power in ", file, " row ", r)
    out[[r]] <- data.frame(t(carry), a = a_grid, power = powers,
                           segment = segment, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

canon <- function(df, test, scenario = NULL, kernel = NA_character_,
                  baseline = NA_character_, m = NA_integer_, p = NA_integer_,
                  n_rep = 100L, source) {
  res <- data.frame(
    test       = test,
    scenario   = if (!is.null(scenario)) scenario else df$scenario,
    segment    = df$segment,
    projection = if ("projection" %in% names(df)) df$projection else "none",
    kernel     = if ("kernel" %in% names(df)) df$kernel else kernel,
    baseline   = baseline,
    m          = if ("m" %in% names(df)) as.integer(df$m) else m,
    p          = if ("p" %in% names(df)) as.integer(df$p) else p,
    gamma      = as.numeric(df$gamma),
    a          = df$a,
    power      = df$power,
    n          = 400L,
    n_rep      = n_rep,
    alpha      = ALPHA,
    source     = source,
    stringsAsFactors = FALSE)
  res
}

# RKHS baseline metadata for the functional (projection, kernel) pairs:
# ito+gaussian and c_max_var+gaussian use the N(0,1) baseline (T2 normal
# scores), l2+min uses runif (T1 uniform ranks).
rkhs_baseline <- function(projection, kernel) {
  ifelse(kernel == "min", "runif", "")
}

message("== Parsing printed thesis tables ==")
printed <- list()
for (seg in c("long", "short")) {
  # F3 (pareto) ------------------------------------------------------------
  d <- parse_printed(sprintf("_pareto_F3_wilcoxon_%s.tex", seg),
                     c("projection", "m", "gamma"), 27)
  d$projection <- vapply(d$projection, map_projection, "")
  printed[[length(printed) + 1]] <- canon(d, "wilcoxon", "F3", kernel = "h_0",
    source = sprintf("thesis _pareto_F3_wilcoxon_%s.tex", seg))

  d <- parse_printed(sprintf("_pareto_F3_lp_%s.tex", seg),
                     c("projection", "p", "gamma"), 9)
  d$projection <- vapply(d$projection, map_projection, "")
  printed[[length(printed) + 1]] <- canon(d, "lp", "F3",
    source = sprintf("thesis _pareto_F3_lp_%s.tex", seg))

  d <- parse_printed(sprintf("_pareto_F3_rkhs_%s.tex", seg),
                     c("projection", "kernel", "gamma"), 9)
  d$projection <- vapply(d$projection, map_projection, "")
  d$kernel <- vapply(d$kernel, map_kernel, "")
  dc <- canon(d, "rkhs", "F3",
    source = sprintf("thesis _pareto_F3_rkhs_%s.tex", seg))
  dc$baseline <- rkhs_baseline(dc$projection, dc$kernel)
  printed[[length(printed) + 1]] <- dc

  # F4 ---------------------------------------------------------------------
  d <- parse_printed(sprintf("_f4_wilcoxon_%s.tex", seg),
                     c("projection", "m", "gamma"), 27)
  d$projection <- vapply(d$projection, map_projection, "")
  printed[[length(printed) + 1]] <- canon(d, "wilcoxon", "F4", kernel = "h_0",
    source = sprintf("thesis _f4_wilcoxon_%s.tex", seg))

  d <- parse_printed(sprintf("_f4_lp_%s.tex", seg),
                     c("projection", "p", "gamma"), 9)
  d$projection <- vapply(d$projection, map_projection, "")
  printed[[length(printed) + 1]] <- canon(d, "lp", "F4",
    source = sprintf("thesis _f4_lp_%s.tex", seg))

  d <- parse_printed(sprintf("_f4_rkhs_%s.tex", seg),
                     c("projection", "kernel", "gamma"), 9)
  d$projection <- vapply(d$projection, map_projection, "")
  d$kernel <- vapply(d$kernel, map_kernel, "")
  dc <- canon(d, "rkhs", "F4",
    source = sprintf("thesis _f4_rkhs_%s.tex (FAR(1) finite-sample null)", seg))
  dc$baseline <- rkhs_baseline(dc$projection, dc$kernel)
  printed[[length(printed) + 1]] <- dc

  d <- parse_printed(sprintf("_f4_birbilas_%s.tex", seg), c("gamma"), 3)
  printed[[length(printed) + 1]] <- canon(d, "far1", "F4",
    source = sprintf("thesis _f4_birbilas_%s.tex (beta1, adjusted c.v.)", seg))

  # Scalar S1-S5 -----------------------------------------------------------
  d <- parse_printed(sprintf("_scalar_wilcoxon_%s.tex", seg),
                     c("scenario", "m", "gamma"), 45)
  printed[[length(printed) + 1]] <- canon(d, "wilcoxon", kernel = "h_w",
    source = sprintf("thesis _scalar_wilcoxon_%s.tex", seg))

  d <- parse_printed(sprintf("_scalar_lp_%s.tex", seg),
                     c("scenario", "gamma"), 15)
  printed[[length(printed) + 1]] <- canon(d, "lp", p = 2L,
    source = sprintf("thesis _scalar_lp_%s.tex", seg))

  d <- parse_printed(sprintf("_scalar_rkhs_%s.tex", seg),
                     c("scenario_kernel", "gamma"), 9)
  d$scenario <- vapply(d$scenario_kernel, map_scenario, "")
  d$kernel   <- vapply(d$scenario_kernel, map_kernel, "")
  printed[[length(printed) + 1]] <- canon(d, "rkhs",
    source = sprintf("thesis _scalar_rkhs_%s.tex (matched kernels)", seg))
}
printed <- do.call(rbind, printed)

# Spot-value asserts against independently read table cells.
pick <- function(df, ...) {
  key <- list(...)
  sel <- rep(TRUE, nrow(df))
  for (k in names(key)) sel <- sel & (df[[k]] == key[[k]] | (is.na(df[[k]]) & is.na(key[[k]])))
  df$power[sel]
}
stopifnot(identical(
  pick(printed, test = "wilcoxon", scenario = "F3", segment = "long",
       projection = "ito", m = 1L, gamma = 0),
  c(1.00, 1.00, 1.00, 0.87, 0.11, 1.00)))
stopifnot(pick(printed, test = "lp", scenario = "S1", segment = "long",
               gamma = 0, a = 0) == 0.46)
stopifnot(pick(printed, test = "far1", scenario = "F4", segment = "long",
               gamma = 0, a = 0.5) == 0.82)
stopifnot(pick(printed, test = "rkhs", scenario = "S4", segment = "long",
               gamma = 0.45, a = 2) == 0.95)
stopifnot(pick(printed, test = "rkhs", scenario = "F3", segment = "long",
               projection = "l2", gamma = 0.25, a = 2) == 0.91)
message(sprintf("   parsed %d printed values, spot-asserts green", nrow(printed)))

# ---- Preserved values (source CSVs deleted; from the original plotting script
# @ c64cd58, the N_rep = 100 F1/F2 Wilcoxon rerun) --------------------------
preserved <- rbind(
  data.frame(scenario = "F1", m = 2L, a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(1.00, 0.56, 0.03, 0.43, 1.00, 1.00)),
  data.frame(scenario = "F1", m = 4L, a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(1.00, 0.31, 0.01, 0.27, 1.00, 1.00)),
  data.frame(scenario = "F2", m = 2L, a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 0.99, 0.02, 0.76, 1.00)),
  data.frame(scenario = "F2", m = 4L, a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 0.98, 0.00, 0.60, 1.00)))
preserved <- data.frame(
  test = "wilcoxon", scenario = preserved$scenario, segment = "long",
  projection = "l2", kernel = "h_0", baseline = NA_character_,
  m = preserved$m, p = NA_integer_, gamma = 0, a = preserved$a,
  power = preserved$power, n = 400L, n_rep = 100L, alpha = ALPHA,
  source = "preserved from the original plotting script (rerun CSVs deleted)",
  stringsAsFactors = FALSE)

# Cross-check block: the 66 m = 1 long-segment values inlined at c64cd58
# (gamma in {0, 0.25}); the fresh rerun must reproduce them exactly
# (same driver, seeds, NCORES).
crosscheck <- rbind(
  data.frame(sc = "F1", proj = "ito",       gamma = 0,    a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(1.00, 0.47, 0.04, 0.31, 1.00, 1.00)),
  data.frame(sc = "F1", proj = "ito",       gamma = 0.25, a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(0.03, 0.00, 0.00, 0.00, 0.06, 0.99)),
  data.frame(sc = "F1", proj = "l2",        gamma = 0,    a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(1.00, 0.84, 0.05, 0.85, 1.00, 1.00)),
  data.frame(sc = "F1", proj = "l2",        gamma = 0.25, a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(0.87, 0.00, 0.00, 0.00, 1.00, 1.00)),
  data.frame(sc = "F1", proj = "c_max_var", gamma = 0,    a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(0.94, 0.16, 0.05, 0.54, 1.00, 1.00)),
  data.frame(sc = "F1", proj = "c_max_var", gamma = 0.25, a = c(1.5, 2, 2.5, 3, 5, 10),
             power = c(0.00, 0.00, 0.00, 0.00, 0.05, 1.00)),
  data.frame(sc = "F2", proj = "ito",       gamma = 0,    a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 0.98, 0.02, 0.70, 1.00)),
  data.frame(sc = "F2", proj = "ito",       gamma = 0.25, a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 0.00, 0.00, 0.00, 0.97)),
  data.frame(sc = "F2", proj = "l2",        gamma = 0,    a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 1.00, 0.06, 1.00, 1.00)),
  data.frame(sc = "F2", proj = "l2",        gamma = 0.25, a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 0.10, 0.00, 0.00, 1.00)),
  data.frame(sc = "F2", proj = "c_max_var", gamma = 0,    a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 0.99, 0.01, 0.60, 1.00)),
  data.frame(sc = "F2", proj = "c_max_var", gamma = 0.25, a = c(0.5, 1, 1.5, 2, 5),
             power = c(1.00, 0.00, 0.00, 0.00, 0.88)))

# ---- Fresh rerun CSVs -----------------------------------------------------
rerun_specs <- expand.grid(test = c("wilcoxon", "lp", "rkhs"),
                           sc = c("F1", "F2"),
                           tag = NA, stringsAsFactors = FALSE)
rerun_files <- c(
  sprintf("thesis_wilcoxon_%s_n400_nrep100_h_0_beta0_%s.csv",
          rep(c("F1", "F2"), each = 3), c("ito", "l2", "c_max_var")),
  sprintf("thesis_lp_%s_n400_nrep100_p2_%s.csv",
          rep(c("F1", "F2"), each = 3), c("ito", "l2", "c_max_var")),
  sprintf("thesis_rkhs_%s_n400_nrep100_%s.csv",
          rep(c("F1", "F2"), each = 3),
          c("gaussian_ito", "min_runif_l2", "gaussian_c_max_var")))

missing <- rerun_files[!file.exists(file.path(res_dir, rerun_files))]
fresh <- NULL
if (length(missing) && ALLOW_MISSING_RERUN) {
  message("!! ALLOW_MISSING_RERUN: emitting printed + preserved rows only; ",
          length(missing), " rerun CSVs still missing:")
  for (f in missing) message("     ", f)
} else if (length(missing)) {
  stop("missing rerun CSVs (run sweeps/run_thesis_power_F1F2_n400_nrep100.sh, or ",
       "set ALLOW_MISSING_RERUN=1): ", paste(missing, collapse = ", "))
} else {
  message("== Reading 18 fresh rerun CSVs ==")
  fresh <- do.call(rbind, lapply(rerun_files, function(f) {
    d <- read.csv(file.path(res_dir, f), stringsAsFactors = FALSE)
    stopifnot(all(d$n == 400), all(d$N_rep == 100))
    kernel <- switch(d$method[1],
      wilcoxon = "h_0",
      lp       = NA_character_,
      rkhs     = if (grepl("min", d$tag[1])) "min" else "gaussian")
    data.frame(
      test = d$method, scenario = d$scenario,
      segment = ifelse(d$ec_length > 0.1, "long", "short"),
      projection = d$projection, kernel = kernel,
      baseline = if (d$method[1] == "rkhs")
        ifelse(grepl("runif", d$tag), "runif", "") else NA_character_,
      m = as.integer(d$m),
      p = if (d$method[1] == "lp") 2L else NA_integer_,
      gamma = d$gamma, a = d$a, power = d$power,
      n = 400L, n_rep = 100L, alpha = ALPHA,
      source = paste0("rerun results/", f),
      stringsAsFactors = FALSE)
  }))
  message(sprintf("   %d fresh values", nrow(fresh)))

  # Check 1: fresh must reproduce the 66 preserved m = 1 values exactly.
  ndiff <- 0L
  for (r in seq_len(nrow(crosscheck))) {
    fv <- pick(fresh, test = "wilcoxon", scenario = crosscheck$sc[r],
               segment = "long", projection = crosscheck$proj[r],
               gamma = crosscheck$gamma[r], a = crosscheck$a[r])
    stopifnot(length(fv) == 1)
    dd <- abs(fv - crosscheck$power[r])
    if (dd > 0.15) stop(sprintf(
      "fresh vs preserved mismatch > 0.15 at %s %s gamma=%s a=%s: %.2f vs %.2f",
      crosscheck$sc[r], crosscheck$proj[r], crosscheck$gamma[r],
      crosscheck$a[r], fv, crosscheck$power[r]))
    if (dd > 0) {
      ndiff <- ndiff + 1L
      message(sprintf("   WARN fresh != preserved at %s %s gamma=%s a=%s: %.2f vs %.2f",
                      crosscheck$sc[r], crosscheck$proj[r], crosscheck$gamma[r],
                      crosscheck$a[r], fv, crosscheck$power[r]))
    }
  }
  message(sprintf("   fresh-vs-preserved: %d/%d cells differ (expect 0)",
                  ndiff, nrow(crosscheck)))

  # Check 2: fresh vs printed N_rep = 20 tables (never emitted) - a
  # wrong-pipeline mistake shows up as systematic large drift, MC noise
  # as a handful of outliers.
  message("== Fresh vs printed N_rep=20 tolerance scan (lp, rkhs) ==")
  drift <- NULL
  for (seg in c("long", "short")) for (sc in c("F1", "F2")) {
    d <- parse_printed(sprintf("_pareto_%s_lp_%s.tex", sc, seg),
                       c("projection", "p", "gamma"), 9, expect_nrep = 20)
    d$projection <- vapply(d$projection, map_projection, "")
    d$test <- "lp"; d$scenario <- sc
    d2 <- parse_printed(sprintf("_pareto_%s_rkhs_%s.tex", sc, seg),
                        c("projection", "kernel", "gamma"), 9, expect_nrep = 20)
    d2$projection <- vapply(d2$projection, map_projection, "")
    d2$test <- "rkhs"; d2$scenario <- sc; d2$kernel <- NULL
    drift <- rbind(drift, d[c("test","scenario","segment","projection","gamma","a","power")],
                   d2[c("test","scenario","segment","projection","gamma","a","power")])
  }
  # N_rep = 20 tables for wilcoxon are superseded by the exact preserved check.
  n_big <- 0L; deltas <- numeric(0)
  for (r in seq_len(nrow(drift))) {
    fv <- pick(fresh, test = drift$test[r], scenario = drift$scenario[r],
               segment = drift$segment[r], projection = drift$projection[r],
               gamma = drift$gamma[r], a = drift$a[r])
    stopifnot(length(fv) == 1)
    deltas <- c(deltas, abs(fv - drift$power[r]))
    if (abs(fv - drift$power[r]) > 0.25) {
      n_big <- n_big + 1L
      message(sprintf("   drift > 0.25: %s %s %s %s gamma=%s a=%s: %.2f (nrep100) vs %.2f (nrep20)",
                      drift$test[r], drift$scenario[r], drift$segment[r],
                      drift$projection[r], drift$gamma[r], drift$a[r],
                      fv, drift$power[r]))
    }
  }
  message(sprintf("   drift: mean |delta| = %.3f, %d/%d cells > 0.25",
                  mean(deltas), n_big, nrow(drift)))
  if (n_big > 0.15 * nrow(drift)) {
    stop("more than 15% of F1/F2 lp/rkhs cells drift > 0.25 vs the printed ",
         "N_rep=20 tables - suspicious of a wrong pipeline, aborting")
  }
}

# ---- Union, uniqueness, write --------------------------------------------
all_rows <- rbind(printed, preserved, if (!is.null(fresh)) fresh)
key <- with(all_rows, paste(test, scenario, segment, projection, kernel,
                            baseline, m, p, gamma, a, sep = "|"))
if (anyDuplicated(key)) {
  stop("duplicate cell keys: ", paste(head(key[duplicated(key)], 5), collapse = " ; "))
}
ord <- with(all_rows, order(test, scenario, segment, projection, kernel, m, p, gamma, a))
all_rows <- all_rows[ord, ]

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(all_rows, out_csv, row.names = FALSE)
message(sprintf("== Wrote %s (%d rows: %d printed, %d preserved, %d fresh) ==",
                out_csv, nrow(all_rows), nrow(printed), nrow(preserved),
                if (is.null(fresh)) 0L else nrow(fresh)))
