rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(ggvenn)
library(dplyr)
library(patchwork)

# ==========================================
# PATHS
# ==========================================
path    <- "~/Documents/projects/ISF/ISF_fede"
deg_dir <- file.path(path, "code", "DEGs", "res")
dirOut  <- file.path(path, "code", "supplementary_analysis", "output")
dir.create(dirOut, showWarnings = FALSE)

# ==========================================
# LOAD DEGs (already filtered FDR <= 0.05)
# ==========================================
deg_IF    <- read.delim(file.path(deg_dir, "IF",    "DEG_res", "DEG.txt"), sep = "\t")
deg_SERUM <- read.delim(file.path(deg_dir, "SERUM", "DEG_res", "DEG.txt"), sep = "\t")

shared_proteins <- intersect(deg_IF$protein, deg_SERUM$protein)
if_unique       <- deg_IF[!(deg_IF$protein %in% deg_SERUM$protein), ]
serum_unique    <- deg_SERUM[!(deg_SERUM$protein %in% deg_IF$protein), ]

message("IF total DEGs: ",    nrow(deg_IF),    " | IF-unique: ",    nrow(if_unique))
message("SERUM total DEGs: ", nrow(deg_SERUM), " | SERUM-unique: ", nrow(serum_unique))
message("Shared: ", length(shared_proteins), " (", paste(shared_proteins, collapse = ", "), ")")

# ==========================================
# FIGURE 1A — Venn diagram
# ==========================================
venn_list <- list(
  IF    = deg_IF$protein,
  SERUM = deg_SERUM$protein
)

p_venn <- ggvenn(
  venn_list,
  fill_color      = c("#2196F3", "#FF5722"),
  stroke_color    = c("#1565C0", "#BF360C"),
  set_name_color  = c("#1565C0", "#BF360C"),
  show_percentage = FALSE,
  text_size       = 5
) +
  labs(title    = "DEGs: T2D vs Control",
       subtitle = "FDR <= 0.05") +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

# ==========================================
# COLOUR PALETTE — direction
# ==========================================
dir_colors <- c(
  "up_in_T2D"   = "#E53935",   # red
  "down_in_T2D" = "#1E88E5"    # blue
)

dir_labels <- c(
  "up_in_T2D"   = "Up in T2D",
  "down_in_T2D" = "Down in T2D"
)

# ==========================================
# FIGURE 1B — IF-unique lollipop
# ==========================================
df_if <- if_unique
df_if$protein <- factor(df_if$protein,
                         levels = df_if$protein[order(df_if$log2FC)])

p_if <- ggplot(df_if, aes(x = log2FC, y = protein, color = direction)) +
  geom_segment(aes(xend = 0, yend = protein), linewidth = 0.9) +
  geom_point(aes(size = -log10(FDR))) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_color_manual(values = dir_colors, labels = dir_labels,
                     name = "Direction") +
  scale_size_continuous(range = c(3, 8), name = "-log10(FDR)") +
  labs(title    = "IF-unique DEGs",
       subtitle = paste0(nrow(df_if), " proteins | not significant in SERUM"),
       x = "log2 Fold Change (T2D vs Control)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", color = "#1565C0"),
        legend.position    = "right")

# ==========================================
# FIGURE 1C — SERUM-unique lollipop
# ==========================================
df_se <- serum_unique
df_se$protein <- factor(df_se$protein,
                         levels = df_se$protein[order(df_se$log2FC)])

p_se <- ggplot(df_se, aes(x = log2FC, y = protein, color = direction)) +
  geom_segment(aes(xend = 0, yend = protein), linewidth = 0.9) +
  geom_point(aes(size = -log10(FDR))) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_color_manual(values = dir_colors, labels = dir_labels,
                     name = "Direction") +
  scale_size_continuous(range = c(3, 8), name = "-log10(FDR)") +
  labs(title    = "SERUM-unique DEGs",
       subtitle = paste0(nrow(df_se), " proteins | not significant in IF"),
       x = "log2 Fold Change (T2D vs Control)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", color = "#BF360C"),
        legend.position    = "right")

# ==========================================
# SAVE
# ==========================================
ggsave(file.path(dirOut, "fig1a_venn_DEGs.pdf"), p_venn, width = 5, height = 5)

p_lolli <- p_if + p_se + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")
ggsave(file.path(dirOut, "fig1b_unique_DEGs_lollipop.pdf"), p_lolli, width = 14, height = 7)

p_full <- (p_venn / plot_spacer()) + (p_if + p_se) +
  plot_layout(ncol = 2, widths = c(1, 2.5)) +
  plot_annotation(
    title    = "DEG Analysis: T2D vs Control - IF vs SERUM",
    subtitle = "IF-unique DEGs: vascular/ECM proteins | SERUM-unique DEGs: cytokines/immune mediators",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )
ggsave(file.path(dirOut, "fig1_DEG_overview.pdf"), p_full, width = 18, height = 9)

message("\nFigure 1 saved to: ", dirOut)
