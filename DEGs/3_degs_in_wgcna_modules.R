rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)

# ==========================================
# PATHS
# ==========================================
path      <- "~/Documents/projects/ISF/ISF_fede"
dirRes    <- file.path(path, "code", "DEGs", "res")
dirCommon <- file.path(path, "code", "DEGs", "common")
dirOut    <- dirCommon

wgcna_base <- file.path(path, "code", "WGCNA", "code", "project", "ISF", "dataset")

path_geneInfo_IF    <- file.path(wgcna_base, "IF",    "Results", "txtFile", "geneInfo.txt")
path_geneInfo_SERUM <- file.path(wgcna_base, "SERUM", "Results", "txtFile", "geneInfo.txt")

fdr_thr <- 0.05

# ==========================================
# Carica tutti i DEG (FDR <= 0.05)
# ==========================================
deg_if <- read.table(file.path(dirRes, "IF", "DEG.txt"),
                     header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
deg_serum <- read.table(file.path(dirRes, "SERUM", "DEG.txt"),
                        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

sig_if    <- deg_if[deg_if$FDR <= fdr_thr, ]
sig_serum <- deg_serum[deg_serum$FDR <= fdr_thr, ]

common_proteins <- intersect(sig_if$protein, sig_serum$protein)

# ==========================================
# Carica geneInfo WGCNA
# ==========================================
load_geneInfo <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", quote = "",
                   check.names = FALSE, row.names = 1)
  rownames(df) <- gsub('"', '', trimws(rownames(df)))
  df
}

gi_if    <- load_geneInfo(path_geneInfo_IF)
gi_serum <- load_geneInfo(path_geneInfo_SERUM)

# ==========================================
# Tabella completa: tutti i DEG IF
# ==========================================
all_if <- data.frame(
  protein   = sig_if$protein,
  log2FC    = sig_if$log2FC,
  FDR       = sig_if$FDR,
  direction = sig_if$direction,
  module    = gi_if[sig_if$protein, "moduleColor"],
  common    = sig_if$protein %in% common_proteins,
  stringsAsFactors = FALSE
)
all_if <- all_if[order(all_if$FDR), ]

# ==========================================
# Tabella completa: tutti i DEG SERUM
# ==========================================
all_serum <- data.frame(
  protein   = sig_serum$protein,
  log2FC    = sig_serum$log2FC,
  FDR       = sig_serum$FDR,
  direction = sig_serum$direction,
  module    = gi_serum[sig_serum$protein, "moduleColor"],
  common    = sig_serum$protein %in% common_proteins,
  stringsAsFactors = FALSE
)
all_serum <- all_serum[order(all_serum$FDR), ]

# ==========================================
# Salva tabelle
# ==========================================
write.table(all_if,
            file = file.path(dirOut, "all_DEGs_IF_wgcna_modules.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(all_serum,
            file = file.path(dirOut, "all_DEGs_SERUM_wgcna_modules.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\n=== DEG IF ===")
print(all_if[, c("protein", "direction", "log2FC", "FDR", "module", "common")])

message("\n=== DEG SERUM ===")
print(all_serum[, c("protein", "direction", "log2FC", "FDR", "module", "common")])

# ==========================================
# Barplot per modulo con direzione (IF)
# ==========================================
make_barplot <- function(df, fluid) {
  df$module[is.na(df$module)] <- "not_in_WGCNA"
  p <- ggplot(df, aes(x = module, fill = direction)) +
    geom_bar(position = "stack") +
    scale_fill_manual(values = c("up_in_T2D" = "red", "down_in_T2D" = "steelblue")) +
    labs(title = paste0("DEG ", fluid, " per modulo WGCNA"),
         x = "Modulo", y = "N proteine", fill = "Direzione") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(dirOut, paste0("barplot_modules_", fluid, ".pdf")),
         plot = p, width = 6, height = 4)
}

make_barplot(all_if,    "IF")
make_barplot(all_serum, "SERUM")

# ==========================================
# Lollipop plot: proteina x log2FC, colore = modulo
# ==========================================
make_lollipop <- function(df, fluid) {
  df$module[is.na(df$module)] <- "not_in_WGCNA"

  # moduli presenti -> usa colori WGCNA dove possibile
  module_colors <- c(
    turquoise    = "turquoise3",
    blue         = "royalblue",
    brown        = "saddlebrown",
    yellow       = "gold2",
    green        = "forestgreen",
    grey         = "grey50",
    not_in_WGCNA = "black"
  )
  present_colors <- module_colors[names(module_colors) %in% df$module]

  # ordina per log2FC
  df$protein <- factor(df$protein, levels = df$protein[order(df$log2FC)])

  p <- ggplot(df, aes(x = log2FC, y = protein, color = module)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = log2FC, y = protein, yend = protein),
                 linewidth = 0.8) +
    geom_point(aes(size = -log10(FDR)), shape = 16) +
    scale_color_manual(values = present_colors, name = "WGCNA module") +
    scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
    labs(title  = paste0("DEG ", fluid, ": T2D vs Control"),
         subtitle = "Colore = modulo WGCNA | Dimensione = -log10(FDR)",
         x = "log2FC (T2D vs Control)",
         y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_line(color = "grey92"),
          panel.grid.major.x = element_blank(),
          legend.position    = "right")

  ggsave(file.path(dirOut, paste0("lollipop_DEGs_", fluid, ".pdf")),
         plot = p, width = 7, height = max(4, nrow(df) * 0.35))
  ggsave(file.path(dirOut, paste0("lollipop_DEGs_", fluid, ".png")),
         plot = p, width = 7, height = max(4, nrow(df) * 0.35), dpi = 300)
}

make_lollipop(all_if,    "IF")
make_lollipop(all_serum, "SERUM")

# ==========================================
# Plot combinato IF + SERUM (per confronto diretto su una slide)
# ==========================================
all_if$fluid    <- "IF"
all_serum$fluid <- "SERUM"
combined <- rbind(all_if, all_serum)
combined$module[is.na(combined$module)] <- "not_in_WGCNA"

module_colors_all <- c(
  turquoise    = "turquoise3",
  blue         = "royalblue",
  brown        = "saddlebrown",
  yellow       = "gold2",
  green        = "forestgreen",
  grey         = "grey50",
  not_in_WGCNA = "black"
)
present_all <- module_colors_all[names(module_colors_all) %in% combined$module]

# ordine proteine per log2FC medio
prot_order <- combined |>
  (\(d) tapply(d$log2FC, d$protein, mean))() |>
  sort() |>
  names()
combined$protein <- factor(combined$protein, levels = prot_order)

p_combined <- ggplot(combined, aes(x = log2FC, y = protein, color = module)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = log2FC, y = protein, yend = protein),
               linewidth = 0.8) +
  geom_point(aes(size = -log10(FDR)), shape = 16) +
  scale_color_manual(values = present_all, name = "WGCNA module") +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
  facet_wrap(~ fluid, ncol = 2, scales = "free_y") +
  labs(title    = "DEG T2D vs Control: IF e SERUM",
       subtitle = "Colore = modulo WGCNA | Dimensione = -log10(FDR)",
       x = "log2FC (T2D vs Control)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_line(color = "grey92"),
        panel.grid.major.x = element_blank(),
        legend.position    = "right",
        strip.text         = element_text(face = "bold", size = 13))

n_prot <- length(unique(combined$protein))
ggsave(file.path(dirOut, "lollipop_DEGs_combined.pdf"),
       plot = p_combined, width = 12, height = max(5, n_prot * 0.3))
ggsave(file.path(dirOut, "lollipop_DEGs_combined.png"),
       plot = p_combined, width = 12, height = max(5, n_prot * 0.3), dpi = 300)

message("\nOutput in: ", dirOut)
