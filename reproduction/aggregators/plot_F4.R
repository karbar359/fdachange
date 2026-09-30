#!/usr/bin/env Rscript
# Power-curve plots for the F4 cross-method comparison. Reads
# `results/F4_combined.csv` (produced by aggregate_F4.R) and writes PNGs
# under `results/F4_plots/`.
#
# Three plots:
#   1. F4_power_cross_method.png   one line per (method, projection)
#                                  pair, faceted by gamma x ec_length.
#                                  Headline for the supervisor.
#   2. F4_power_wilcoxon_m.png     Wilcoxon M in {1, 2, 4} per
#                                  projection. Shows the m-split's
#                                  effect on power at gamma = 0.
#   3. F4_size_at_a0.png           size at alpha = 0.05 across cells
#                                  and gammas; horizontal reference at
#                                  the nominal level. Surfaces the
#                                  Lp/Wilcoxon over-rejection on FAR(1)
#                                  data at a glance.
#
# Usage:
#   Rscript reproduction/plot_F4.R

suppressPackageStartupMessages(library(ggplot2))

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
res_dir <- file.path(root, "results")
out_dir <- file.path(res_dir, "F4_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(file.path(res_dir, "F4_combined.csv"), stringsAsFactors = FALSE)
df$ec_label <- sprintf("l/n = %.2f", df$ec_length)
df$gamma_label <- sprintf("gamma = %.2f", df$gamma)

ALPHA <- 0.05
A_GRID <- sort(unique(df$a))
plot_theme <- theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())

# ---- 1. Cross-method headline (canonical pair per method) ---------------
# One curve per method, picking the projection that gives the cleanest
# size + power tradeoff so the supervisor can read it.
canon <- rbind(
  df[df$method == "far1", ],
  df[df$method == "lp"       & df$projection == "c_max_var", ],
  df[df$method == "wilcoxon" & df$projection == "c_max_var" & df$m == 1L, ],
  df[df$method == "rkhs"     & grepl("^min_rnorm", df$tag), ])
canon$method_label <- factor(canon$method,
  levels = c("far1", "rkhs", "lp", "wilcoxon"),
  labels = c("Birbilas E (adjusted CV)", "RKHS (min, l2)",
             "Lp (p=2, c_max_var)", "Wilcoxon (h_0, c_max_var, M=1)"))

p1 <- ggplot(canon,
             aes(x = a, y = power,
                 colour = method_label, group = method_label)) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = A_GRID, trans = "sqrt") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_grid(rows = vars(ec_label), cols = vars(gamma_label)) +
  labs(title = "F4 cross-method power curves (canonical pair per method)",
       subtitle = "n = 200, NREP = 100. alpha = 0.05 (dotted). x-axis sqrt scale.",
       x = "jump magnitude a", y = "rejection rate",
       colour = NULL) +
  plot_theme
ggsave(file.path(out_dir, "F4_power_cross_method.png"), p1,
       width = 11, height = 7, dpi = 150, bg = "white")

# ---- 1b. All projections per method (auxiliary, fuller picture) --------
df_m1 <- df[df$m == 1L, ]
df_m1$method_proj <- paste0(df_m1$method, " : ",
                            ifelse(df_m1$projection == "none",
                                   "none", df_m1$projection))
p1b <- ggplot(df_m1,
              aes(x = a, y = power,
                  colour = projection, linetype = projection,
                  group = method_proj)) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = A_GRID, trans = "sqrt") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_grid(rows = vars(method, ec_label), cols = vars(gamma_label)) +
  labs(title = "F4 power curves: all (method, projection) cells (M = 1)",
       subtitle = "n = 200, NREP = 100. alpha = 0.05 (dotted). x-axis sqrt scale.",
       x = "jump magnitude a", y = "rejection rate") +
  plot_theme
ggsave(file.path(out_dir, "F4_power_all_projections.png"), p1b,
       width = 12, height = 12, dpi = 150, bg = "white")

# ---- 2. Wilcoxon M-sweep ------------------------------------------------
df_w <- df[df$method == "wilcoxon", ]
df_w$M_label <- factor(sprintf("M = %d", df_w$m),
                       levels = paste0("M = ", c(1, 2, 4)))

p2 <- ggplot(df_w,
             aes(x = a, y = power,
                 colour = M_label, group = M_label)) +
  geom_hline(yintercept = ALPHA, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  scale_x_continuous(breaks = A_GRID, trans = "sqrt") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_grid(rows = vars(projection, ec_label),
             cols = vars(gamma_label)) +
  labs(title = "F4 Wilcoxon: power vs M (T^(d) outer-max splits)",
       subtitle = "n = 200, kernel h_0, NREP = 100. alpha = 0.05 (dotted).",
       x = "jump magnitude a", y = "rejection rate",
       colour = "splits") +
  plot_theme
ggsave(file.path(out_dir, "F4_power_wilcoxon_m.png"), p2,
       width = 12, height = 11, dpi = 150, bg = "white")

# ---- 3. Size at a = 0 ---------------------------------------------------
df_size <- df[df$a == 0, ]
df_size$cell <- paste(df_size$method, df_size$tag, sep = ":")
df_size <- df_size[order(df_size$method, df_size$cell), ]
df_size$cell <- factor(df_size$cell, levels = unique(df_size$cell))

p3 <- ggplot(df_size,
             aes(x = cell, y = power,
                 colour = factor(gamma), shape = ec_label)) +
  geom_hline(yintercept = ALPHA, linetype = "dashed", colour = "grey40") +
  geom_point(size = 3.0, position = position_dodge(width = 0.45)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "F4 empirical size at a = 0",
       subtitle = "n = 200, NREP = 100. alpha = 0.05 dashed. Lp / Wilcoxon-at-gamma=0 over-reject (no FAR(1) deweighting).",
       x = NULL, y = "rejection rate at a = 0",
       colour = "gamma", shape = "ec length") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "F4_size_at_a0.png"), p3,
       width = 13, height = 7, dpi = 150, bg = "white")

cat("Wrote 4 plots to ", out_dir, "/\n", sep = "")
