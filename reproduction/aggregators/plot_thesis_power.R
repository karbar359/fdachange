#!/usr/bin/env Rscript
# Power-visualization PNGs for the thesis power figures.
#
# One comparable visual system over the full grid
#   {Wilcoxon, Lp, RKHS} x {Ito, L2 norm, C max-var} x {F1-F4}   (functional)
#   {Wilcoxon, Lp, RKHS} x {S1-S5}                               (scalar)
# plus the F4 cross-method figure (incl. Birbilas) and the m-effect backup.
# Encoding everywhere: lines/points = gamma (blue circle 0, aqua triangle
# 0.25, yellow square 0.45), panels = structure, dotted line = alpha.
#
# Reads ONLY the committed tidy table data/thesis_power_values.csv,
# built by code/assemble_thesis_power.R (provenance in its `source` column).
#
# Usage:
#   Rscript reproduction/plot_thesis_power.R
#   ALLOW_PARTIAL_GRID=1 Rscript reproduction/plot_thesis_power.R   # layout drafts
#                                  # while the F1/F2 rerun is still running

suppressPackageStartupMessages(library(ggplot2))

# ---- Constants (edit here) ---------------------------------------------
ALPHA      <- 0.05     # nominal level (dotted reference line)
DPI        <- 170
PNG_W      <- 10.0     # inches, all figures
H_FUNC     <- 6.0      # functional 3x4 grids
H_SCALAR   <- 5.8      # scalar 3x5 grid
H_F4       <- 4.8      # F4 cross-method 2x4 grid
H_MEFFECT  <- 4.0      # m-effect 1x3 strip
BG         <- "white"
ALLOW_PARTIAL_GRID <- identical(Sys.getenv("ALLOW_PARTIAL_GRID"), "1")

# Validated reference categorical palette (slots 1-3, light mode):
COL1 <- "#2a78d6"  # blue
COL2 <- "#1baf7a"  # aqua
COL3 <- "#eda100"  # yellow
INK_PRIMARY   <- "#0b0b0b"
INK_SECONDARY <- "#52514e"
INK_MUTED     <- "#898781"
GRIDLINE      <- "#e1e0d9"

# ---- Paths -------------------------------------------------------------
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- sub(key, "", args[grepl(key, args)])
  if (length(hit)) return(normalizePath(hit[1]))
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) return(normalizePath(ofile))
  NULL
})()
root    <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data_csv <- file.path(root, "data", "thesis_power_values.csv")
if (!file.exists(data_csv)) {
  stop("missing ", data_csv,
       " - run `Rscript reproduction/assemble_thesis_power.R` first")
}
dat <- read.csv(data_csv, stringsAsFactors = FALSE)
dat$g <- factor(dat$gamma, levels = c(0, 0.25, 0.45))

# ---- Shared scales and themes -------------------------------------------
gamma_cols   <- c(`0` = COL1, `0.25` = COL2, `0.45` = COL3)
gamma_shapes <- c(`0` = 16, `0.25` = 17, `0.45` = 15)
gamma_labels <- c(expression(gamma == 0), expression(gamma == 0.25),
                  expression(gamma == 0.45))
scale_gamma <- list(
  scale_colour_manual(values = gamma_cols, labels = gamma_labels, name = NULL),
  scale_shape_manual(values = gamma_shapes, labels = gamma_labels, name = NULL),
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)))

# Grid theme for the 8-15 panel figures; base theme kept for the m-effect
# strip so that backup figure stays visually unchanged.
plot_theme_grid <- theme_minimal(base_size = 18) +
  theme(legend.position   = "bottom",
        legend.text       = element_text(size = 15, colour = INK_PRIMARY),
        axis.title        = element_text(size = 15, colour = INK_SECONDARY),
        axis.text         = element_text(size = 11, colour = INK_MUTED),
        strip.text        = element_text(size = 14, colour = INK_PRIMARY,
                                         face = "bold"),
        strip.text.y      = element_text(size = 12, colour = INK_PRIMARY,
                                         face = "bold"),
        plot.caption      = element_text(size = 11, colour = INK_MUTED,
                                         hjust = 1),
        plot.background   = element_rect(fill = BG, colour = NA),
        panel.background  = element_rect(fill = BG, colour = NA),
        panel.grid.major  = element_line(colour = GRIDLINE, linewidth = 0.4),
        panel.grid.minor  = element_blank(),
        panel.spacing     = unit(0.7, "lines"))

plot_theme <- theme_minimal(base_size = 20) +
  theme(legend.position   = "bottom",
        legend.title      = element_text(size = 17, colour = INK_SECONDARY),
        legend.text       = element_text(size = 17, colour = INK_PRIMARY),
        axis.title        = element_text(size = 18, colour = INK_SECONDARY),
        axis.text         = element_text(size = 15, colour = INK_MUTED),
        strip.text        = element_text(size = 17, colour = INK_PRIMARY,
                                         face = "bold"),
        plot.caption      = element_text(size = 13, colour = INK_MUTED),
        plot.background   = element_rect(fill = BG, colour = NA),
        panel.background  = element_rect(fill = BG, colour = NA),
        panel.grid.major  = element_line(colour = GRIDLINE, linewidth = 0.4),
        panel.grid.minor  = element_blank(),
        panel.spacing     = unit(1.0, "lines"))

save_png <- function(p, name, w, h) {
  out <- file.path(out_dir, name)
  ggsave(out, p, width = w, height = h, dpi = DPI, bg = BG)
  cat("Wrote ", out, "\n", sep = "")
}

check_grid <- function(ok, msg) {
  if (ok) return(invisible(TRUE))
  if (ALLOW_PARTIAL_GRID) { message("   PARTIAL: ", msg); return(invisible(FALSE)) }
  stop(msg, " (set ALLOW_PARTIAL_GRID=1 for layout drafts)")
}

F_LEVELS <- paste0("F", 1:4)
F_LABELS <- c(F1 = "F1: N(0,1) vs Pa(a)",
              F2 = "F2: Pa(1.5) vs Pa(a)",
              F3 = "F3: random walk, Pa(a)",
              F4 = "F4: FAR(1) epidemic")
S_LEVELS <- paste0("S", 1:5)
S_LABELS <- c(S1 = "S1: mean shift",
              S2 = "S2: noncentr. chi-sq",
              S3 = "S3: Pareto scale",
              S4 = "S4: Beta shape",
              S5 = "S5: Pareto shape")
PROJ_LEVELS <- c("ito", "l2", "c_max_var")
seg_txt <- c(long = "long segment (l* = n/5)",
             short = "short segment (l* = [log n]+1)")

# =========================================================================
# 1) Per-test functional figures: projections x F1-F4, lines = gamma.
# =========================================================================
functional_fig <- function(test_id, seg, fname, proj_labels, caption_head) {
  d <- dat[dat$test == test_id & dat$scenario %in% F_LEVELS &
           dat$segment == seg & (is.na(dat$m) | dat$m == 1), ]
  check_grid(setequal(unique(d$scenario), F_LEVELS) &&
             setequal(unique(d$projection), PROJ_LEVELS) &&
             setequal(unique(d$gamma), c(0, 0.25, 0.45)),
             paste0(test_id, " ", seg, ": functional grid incomplete"))
  d$sc_lab   <- factor(d$scenario, levels = F_LEVELS, labels = F_LABELS)
  d$proj_lab <- factor(d$projection, levels = PROJ_LEVELS,
                       labels = proj_labels)
  p <- ggplot(d, aes(a, power, colour = g, shape = g)) +
    geom_hline(yintercept = ALPHA, linetype = "dotted", colour = INK_MUTED) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(proj_lab ~ sc_lab, scales = "free_x") +
    scale_gamma +
    labs(x = "change magnitude a", y = "rejection rate") +
    plot_theme_grid
  save_png(p, fname, PNG_W, H_FUNC)
}

proj_labels_plain <- c("Ito projection", "L2 norm", "C max-var")
proj_labels_rkhs  <- c("Ito + Gaussian", "L2 + min", "C max-var + Gaussian")

for (seg in c("long", "short")) {
  functional_fig("wilcoxon", seg,
    sprintf("thesis_wilcoxon_functional_%s.png", seg), proj_labels_plain,
    "Wilcoxon scan, h_0 kernel, m = 1 (F1/F2: N_rep = 100 rerun)")
  functional_fig("lp", seg,
    sprintf("thesis_lp_functional_%s.png", seg), proj_labels_plain,
    "Lp scan, p = 2 (F1/F2: N_rep = 100 rerun)")
  functional_fig("rkhs", seg,
    sprintf("thesis_rkhs_functional_%s.png", seg), proj_labels_rkhs,
    paste0("RKHS scan, rank-calibrated (normal scores / uniform ranks); ",
           "F4 uses the FAR(1) finite-sample null (F1/F2: N_rep = 100 rerun)"))
}

# =========================================================================
# 2) Cross-test scalar figures: tests x S1-S5, lines = gamma.
#    RKHS x {S3, S5} panels stay empty: no Pareto-baseline null cache
#    (README.md, same gap as in the printed thesis).
# =========================================================================
scalar_fig <- function(seg, fname) {
  d <- dat[dat$scenario %in% S_LEVELS & dat$segment == seg &
           ((dat$test == "wilcoxon" & dat$m == 1) |
            (dat$test == "lp") |
            (dat$test == "rkhs")), ]
  check_grid(setequal(unique(d$scenario[d$test != "rkhs"]), S_LEVELS) &&
             setequal(unique(d$scenario[d$test == "rkhs"]), c("S1", "S2", "S4")),
             paste0("scalar ", seg, ": grid incomplete"))
  test_levels <- c("wilcoxon", "lp", "rkhs")
  test_labels <- c("Wilcoxon m=1", "Lp p=2", "RKHS matched")
  d$test_lab <- factor(d$test, levels = test_levels, labels = test_labels)
  d$sc_lab   <- factor(d$scenario, levels = S_LEVELS, labels = S_LABELS)

  kernel_ann <- data.frame(
    test_lab = factor("RKHS matched", levels = test_labels),
    sc_lab   = factor(S_LABELS[c("S1", "S2", "S4")], levels = S_LABELS),
    lab      = c("Gaussian", "min", "2 - max"))
  na_ann <- data.frame(
    test_lab = factor("RKHS matched", levels = test_labels),
    sc_lab   = factor(S_LABELS[c("S3", "S5")], levels = S_LABELS),
    lab      = "no Pareto-baseline\nnull cache")

  p <- ggplot(d, aes(a, power, colour = g, shape = g)) +
    geom_hline(yintercept = ALPHA, linetype = "dotted", colour = INK_MUTED) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    geom_text(data = kernel_ann, aes(x = 0, y = 0.95, label = lab),
              inherit.aes = FALSE, hjust = 0, size = 3.4,
              colour = INK_SECONDARY) +
    geom_text(data = na_ann, aes(x = 2, y = 0.5, label = lab),
              inherit.aes = FALSE, size = 3.4, colour = INK_MUTED,
              lineheight = 0.95) +
    facet_grid(test_lab ~ sc_lab) +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 2, 4, 8)) +
    scale_gamma +
    labs(x = "change magnitude a", y = "rejection rate") +
    plot_theme_grid +
    theme(strip.text.x = element_text(size = 12, colour = INK_PRIMARY,
                                      face = "bold"))
  save_png(p, fname, PNG_W, H_SCALAR)
}

scalar_fig("long",  "thesis_scalar_all_long.png")
scalar_fig("short", "thesis_scalar_all_short.png")

# =========================================================================
# 3) F4 cross-method figure: segment rows x method columns, lines = gamma.
#    Scan tests shown at the L2-norm projection (the thesis-highlighted
#    pairing); Birbilas's E_{n,beta1,gamma} (adjusted c.v.) is projection-free.
# =========================================================================
f4_fig <- function(fname) {
  d <- dat[dat$scenario == "F4" &
           ((dat$test == "wilcoxon" & dat$projection == "l2" & dat$m == 1) |
            (dat$test == "lp"       & dat$projection == "l2") |
            (dat$test == "rkhs"     & dat$projection == "l2") |
            (dat$test == "birbilas")), ]
  check_grid(setequal(unique(d$test),
                      c("wilcoxon", "lp", "rkhs", "birbilas")) &&
             setequal(unique(d$segment), c("long", "short")),
             "F4 cross-method grid incomplete")
  method_levels <- c("wilcoxon", "lp", "rkhs", "birbilas")
  method_labels <- c("Wilcoxon (h_0, m = 1)", "Lp (p = 2)",
                     "RKHS (min)", "Birbilas E (adj. CV)")
  d$m_lab <- factor(d$test, levels = method_levels, labels = method_labels)
  d$seg_lab <- factor(d$segment, levels = c("long", "short"),
                      labels = c("long segment", "short segment"))
  p <- ggplot(d, aes(a, power, colour = g, shape = g)) +
    geom_hline(yintercept = ALPHA, linetype = "dotted", colour = INK_MUTED) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(seg_lab ~ m_lab) +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 2, 4, 8)) +
    scale_gamma +
    labs(x = "change magnitude a", y = "rejection rate") +
    plot_theme_grid
  save_png(p, fname, PNG_W, H_F4)
}

f4_fig("thesis_f4_cross_method.png")

# =========================================================================
# 4) Backup: Wilcoxon splitting parameter m on F1-F3 (L2-norm projection,
#    gamma = 0, long segment). Layout unchanged from the previous deck
#    figure; m = 1 from the fresh rerun, F1/F2 m in {2, 4} preserved,
#    F3 from the printed table.
# =========================================================================
wm <- dat[dat$test == "wilcoxon" & dat$scenario %in% c("F1", "F2", "F3") &
          dat$segment == "long" & dat$projection == "l2" & dat$gamma == 0, ]
check_grid(setequal(unique(wm$m), c(1, 2, 4)) &&
           setequal(unique(wm$scenario), c("F1", "F2", "F3")),
           "m-effect grid incomplete")
wm$sc_lab <- factor(wm$scenario, levels = c("F1", "F2", "F3"),
                    labels = F_LABELS[1:3])
wm$m <- factor(wm$m)

p_wm <- ggplot(wm, aes(a, power, colour = m, shape = m)) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", colour = INK_MUTED) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.0) +
  facet_wrap(~ sc_lab, nrow = 1, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_colour_manual(values = c(`1` = COL1, `2` = COL2, `4` = COL3),
                      name = "splits m") +
  scale_shape_manual(values = c(`1` = 16, `2` = 17, `4` = 15),
                     name = "splits m") +
  labs(x = "Pareto tail index a", y = "rejection rate") +
  plot_theme
save_png(p_wm, "thesis_wilcoxon_m_effect.png", PNG_W, H_MEFFECT)

cat("Done.\n")
