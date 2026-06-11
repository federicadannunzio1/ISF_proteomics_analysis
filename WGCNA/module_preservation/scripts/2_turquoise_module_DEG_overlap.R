rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(ggrepel)

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
base_path <- "~/Documents/projects/ISF/ISF_fede"

path_geneInfo <- file.path(base_path, "code", "WGCNA", "code", "project",
                           "ISF", "dataset", "IF", "Results", "txtFile", "geneInfo.txt")
path_deg_if   <- file.path(base_path, "code", "DEGs", "res", "IF", "DEG.txt")
path_deg_se   <- file.path(base_path, "code", "DEGs", "res", "SERUM", "DEG.txt")

out_dir <- file.path(base_path, "code", "WGCNA", "res", "turquoise_module")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fdr_thr <- 0.05

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────
# All 171 IF proteins with their module assignment and kME
gi <- read.table(path_geneInfo, header = TRUE, sep = "\t", quote = "",
                 check.names = FALSE, row.names = 1)
rownames(gi) <- gsub('"', '', trimws(rownames(gi)))

# All IF DEGs (FDR <= 0.05) with log2FC and direction
deg_if <- read.table(path_deg_if, header = TRUE, sep = "\t", quote = "",
                     stringsAsFactors = FALSE)
deg_if <- deg_if[deg_if$FDR <= fdr_thr, ]

# All SERUM DEGs (to flag proteins also significant in serum)
deg_se <- read.table(path_deg_se, header = TRUE, sep = "\t", quote = "",
                     stringsAsFactors = FALSE)
deg_se <- deg_se[deg_se$FDR <= fdr_thr, ]

# ─────────────────────────────────────────────────────────────────────────────
# Focus on turquoise IF module
# ─────────────────────────────────────────────────────────────────────────────
turq <- data.frame(
  protein       = rownames(gi)[gi$moduleColor == "turquoise"],
  kME_turquoise = gi[gi$moduleColor == "turquoise", "MM.turquoise"],
  stringsAsFactors = FALSE
)

# Annotate with DEG status
turq$is_DEG_IF   <- turq$protein %in% deg_if$protein
turq$is_DEG_SERUM <- turq$protein %in% deg_se$protein

# Merge log2FC and FDR from IF DEGs (NA for non-DEG proteins)
turq <- merge(turq,
              deg_if[, c("protein", "log2FC", "FDR", "direction")],
              by = "protein", all.x = TRUE)

# Direction label for plot
turq$status <- ifelse(!turq$is_DEG_IF, "not DEG",
                      ifelse(turq$direction == "up_in_T2D", "up in T2D", "down in T2D"))
turq$status <- factor(turq$status, levels = c("down in T2D", "not DEG", "up in T2D"))

turq <- turq[order(turq$kME_turquoise, decreasing = TRUE), ]

# ─────────────────────────────────────────────────────────────────────────────
# Summary table
# ─────────────────────────────────────────────────────────────────────────────
message("=== Turquoise module IF: ", nrow(turq), " proteins ===")
message("  DEG in IF:   ", sum(turq$is_DEG_IF))
message("  DEG in SERUM: ", sum(turq$is_DEG_SERUM))
message("")
message("--- DEG proteins in turquoise module ---")
print(turq[turq$is_DEG_IF, c("protein", "log2FC", "FDR", "direction",
                               "kME_turquoise", "is_DEG_SERUM")])

write.table(turq,
            file = file.path(out_dir, "turquoise_module_proteins.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Plot 1: kME vs log2FC (all turquoise proteins; DEGs labeled)
# kME = how strongly each protein belongs to the module (hub = high kME)
# log2FC = only available for DEGs; NA for non-DEGs → plotted at x=0
# ─────────────────────────────────────────────────────────────────────────────
turq$log2FC_plot <- ifelse(is.na(turq$log2FC), 0, turq$log2FC)
turq$label       <- ifelse(turq$is_DEG_IF, turq$protein, "")

p1 <- ggplot(turq, aes(x = log2FC_plot, y = kME_turquoise,
                        color = status, label = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(size = is_DEG_IF), alpha = 0.8) +
  scale_size_manual(values = c("TRUE" = 3.5, "FALSE" = 1.8), guide = "none") +
  scale_color_manual(
    values = c("down in T2D" = "steelblue", "not DEG" = "grey70", "up in T2D" = "red"),
    name   = "DEG status"
  ) +
  geom_text_repel(size = 3, fontface = "bold",
                  box.padding = 0.4, point.padding = 0.2,
                  max.overlaps = Inf, show.legend = FALSE) +
  labs(
    title    = "IF turquoise module: DEG status and module membership",
    subtitle = paste0(nrow(turq), " proteins | kME = module membership (hub = high kME) | ",
                      "non-DEG proteins plotted at log2FC = 0"),
    x = "log2FC (T2D vs Control)",
    y = "kME turquoise (module membership)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8, color = "grey40"))

ggsave(file.path(out_dir, "turquoise_kME_vs_log2FC.pdf"), plot = p1, width = 8, height = 6)
ggsave(file.path(out_dir, "turquoise_kME_vs_log2FC.png"), plot = p1, width = 8, height = 6, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Plot 2: barplot of DEG proteins ordered by log2FC
# ─────────────────────────────────────────────────────────────────────────────
deg_turq <- turq[turq$is_DEG_IF, ]
deg_turq$protein <- factor(deg_turq$protein,
                            levels = deg_turq$protein[order(deg_turq$log2FC)])

p2 <- ggplot(deg_turq, aes(x = log2FC, y = protein, fill = direction)) +
  geom_col() +
  geom_text(aes(label = paste0("FDR=", formatC(FDR, format = "e", digits = 1))),
            hjust = ifelse(deg_turq$log2FC[order(deg_turq$log2FC)] > 0, -0.1, 1.1),
            size = 2.8) +
  scale_fill_manual(values = c("down_in_T2D" = "steelblue", "up_in_T2D" = "red"),
                    labels = c("down in T2D", "up in T2D"),
                    name   = NULL) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  labs(
    title    = "DEG proteins in IF turquoise module",
    subtitle = "Module: not preserved in SERUM (Z=1.46) | r(ME, T2D) = -0.33, p=0.006",
    x = "log2FC (T2D vs Control)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8, color = "grey40"))

ggsave(file.path(out_dir, "turquoise_DEG_barplot.pdf"), plot = p2, width = 7, height = 5)
ggsave(file.path(out_dir, "turquoise_DEG_barplot.png"), plot = p2, width = 7, height = 5, dpi = 300)

message("\nOutput saved in: ", out_dir)
message("Done.")
