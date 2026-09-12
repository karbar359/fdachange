#!/usr/bin/env Rscript
# Interval-agreement figure for the application section.
#
# Draws the detected change segment [l_lo, l_hi] of every test x projection
# combination on the Litgrid data (dissertation Table 7.1, no-outlier run;
# values transcribed to data/litgrid_intervals.csv) as horizontal bars
# on the [0, 1] sample axis:
#   - nine bars at each test's locating gamma (Wilcoxon 0.05, Lp 0, RKHS 0),
#     Ito rows with +-sd whiskers on both endpoints (mean over 10 draws);
#   - two lighter bars for the complementary window the C max-var Lp/RKHS
#     rows flag at gamma = 0.25 - the same break seen from the right.
# A dashed line marks the shared right endpoint ~0.60.
#
# Usage: Rscript reproduction/plot_litgrid_intervals.R

# ---- Inputs (edit these) --------------------------------------------------
IN_CSV  <- NULL                    # NULL -> data/litgrid_intervals.csv
OUT_PNG <- NULL                    # NULL -> output/figures/litgrid_intervals.png
DPI     <- 160
BG      <- "white"
W_IN    <- 10.5
H_IN    <- 5.4
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- sub("--file=", "", args[grepl("^--file=", args)])
  if (length(hit)) return(normalizePath(hit[1]))
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) return(normalizePath(ofile))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))
if (!dir.exists(file.path(root, "data"))) root <- file.path(root, "code")
repo <- normalizePath(file.path(root, ".."))

if (is.null(IN_CSV))  IN_CSV  <- file.path(root, "data", "litgrid_intervals.csv")
if (is.null(OUT_PNG)) OUT_PNG <- file.path(root, "output", "figures",
                                           "litgrid_intervals.png")

INK_PRIMARY   <- "#0b0b0b"
INK_SECONDARY <- "#52514e"
INK_MUTED     <- "#898781"
TEST_COLS <- c(wilcoxon = "#7b003f", lp = "#1f6f8b", rkhs = "#3d7a4f")
TEST_LABS <- c(wilcoxon = "Wilcoxon", lp = "Lp (p = 2)", rkhs = "RKHS")
PROJ_LABS <- c(ito = "Ito", l2 = "L2 norm", c_max_var = "C max-var")

d <- read.csv(IN_CSV)

# Locating rows: the gamma at which each test reports the first window.
loc <- rbind(
  d[d$test == "wilcoxon" & d$gamma == 0.05, ],
  d[d$test == "lp"       & d$gamma == 0,    ],
  d[d$test == "rkhs"     & d$gamma == 0,    ])
loc$kind <- "locating"
loc$glab <- ifelse(loc$test == "wilcoxon", "\u03b3 = 0.05", "\u03b3 = 0")

# Complementary rows: C max-var Lp/RKHS at gamma = 0.25 flag [0.60, 1.0].
comp <- d[d$projection == "c_max_var" & d$gamma == 0.25 & d$test %in% c("lp", "rkhs"), ]
comp$kind <- "complementary"
comp$glab <- "\u03b3 = 0.25"

pd <- rbind(loc, comp)
pd$row_lab <- sprintf("%s · %s (%s)",
                      TEST_LABS[pd$test], PROJ_LABS[pd$projection], pd$glab)
row_order <- c(
  sprintf("%s · %s (%s)", TEST_LABS["wilcoxon"], PROJ_LABS, "\u03b3 = 0.05"),
  sprintf("%s · %s (%s)", TEST_LABS["lp"],       PROJ_LABS, "\u03b3 = 0"),
  sprintf("%s · %s (%s)", TEST_LABS["rkhs"],     PROJ_LABS, "\u03b3 = 0"),
  sprintf("%s · %s (%s)", TEST_LABS[c("lp", "rkhs")], PROJ_LABS["c_max_var"], "\u03b3 = 0.25"))
pd$row_lab <- factor(pd$row_lab, levels = rev(row_order))
pd$alpha_v <- ifelse(pd$kind == "locating", 1, 0.45)

p <- ggplot(pd, aes(y = row_lab, colour = test)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.5,
           fill = NA, colour = NA) +
  geom_vline(xintercept = 0.60, linetype = "dashed", colour = INK_MUTED) +
  geom_segment(aes(x = l_lo, xend = l_hi, yend = row_lab, alpha = alpha_v),
               linewidth = 3.2, lineend = "round") +
  geom_point(aes(x = l_lo, alpha = alpha_v), size = 2.1) +
  geom_point(aes(x = l_hi, alpha = alpha_v), size = 2.1) +
  geom_errorbarh(aes(xmin = l_lo - l_lo_sd, xmax = l_lo + l_lo_sd),
                 height = 0.28, linewidth = 0.5, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = l_hi - l_hi_sd, xmax = l_hi + l_hi_sd),
                 height = 0.28, linewidth = 0.5, na.rm = TRUE) +
  annotate("text", x = 0.615, y = length(row_order) + 0.45,
           label = "shared right endpoint ≈ 0.60", hjust = 0,
           size = 3.8, colour = INK_SECONDARY) +
  scale_colour_manual(values = TEST_COLS, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-0.05, 1.05),
                     breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = "position in the sample (fraction of n)", y = NULL,
       caption = "Detected segments, dissertation Table 7.1 (outlier week removed). Ito bars: mean over 10 draws, whiskers +-1 sd.\nLighter bars: the complementary window flagged at \u03b3 = 0.25.") +
  theme_minimal(base_size = 17) +
  theme(axis.title   = element_text(size = 14, colour = INK_SECONDARY),
        axis.text.y  = element_text(size = 12.5, colour = INK_PRIMARY),
        axis.text.x  = element_text(size = 12, colour = INK_MUTED),
        plot.caption = element_text(size = 9.5, colour = INK_MUTED, hjust = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.background = element_rect(fill = BG, colour = NA))

dir.create(dirname(OUT_PNG), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PNG, p, width = W_IN, height = H_IN, dpi = DPI, bg = BG)
cat("written:", OUT_PNG, "\n")
