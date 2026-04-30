rm(list = ls())

library(ggplot2)
library(reshape2)
library(pheatmap)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir  <- "~/Documents/projects/ISF/ISF_fede"
wgcna_dir <- file.path(base_dir, "code/WGCNA_paola/code/project/ISF/dataset")
out_dir   <- file.path(base_dir, "code/modules_analysis")

# ── Load expression matrices (proteins = rows, samples = cols) ─────────────────
if_mat <- read.table(file.path(wgcna_dir, "IF/matrix/matrix.txt"),
                     sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                     row.names = 1)

serum_mat <- read.table(file.path(wgcna_dir, "SERUM/matrix/matrix.txt"),
                        sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                        row.names = 1)

rownames(if_mat)    <- gsub('"', '', rownames(if_mat))
rownames(serum_mat) <- gsub('"', '', rownames(serum_mat))
colnames(if_mat)    <- gsub('"', '', colnames(if_mat))
colnames(serum_mat) <- gsub('"', '', colnames(serum_mat))

# ── Common proteins ────────────────────────────────────────────────────────────
common_proteins <- c("APOM","ENG","TGFBI","DPP4","MET","DNER",
                     "NCAM1","CHL1","TCN2","PLXNB2","TNXB","FETUB","IL7R")

# Keep only proteins present in both matrices (safety check)
common_proteins <- common_proteins[common_proteins %in% rownames(if_mat) &
                                     common_proteins %in% rownames(serum_mat)]
cat("Proteins retained:", paste(common_proteins, collapse = ", "), "\n")

# ── Extract sub-matrices ───────────────────────────────────────────────────────
if_sub    <- if_mat[common_proteins, , drop = FALSE]
serum_sub <- serum_mat[common_proteins, , drop = FALSE]

# ── Align paired samples (IF X ↔ Serum X) ─────────────────────────────────────
# Extract numeric ID from column names
if_ids    <- as.integer(gsub("^IF\\s+", "", colnames(if_sub)))
serum_ids <- as.integer(gsub("^Serum\\s+", "", colnames(serum_sub)))

shared_ids <- sort(intersect(if_ids, serum_ids))
cat("Paired samples (n):", length(shared_ids), "\n")

if_sub    <- if_sub[, match(shared_ids, if_ids), drop = FALSE]
serum_sub <- serum_sub[, match(shared_ids, serum_ids), drop = FALSE]

# Rename columns to bare sample IDs for legibility
colnames(if_sub)    <- as.character(shared_ids)
colnames(serum_sub) <- as.character(shared_ids)

# ── 1. Side-by-side heatmaps ──────────────────────────────────────────────────
# Scale each matrix row-wise (z-score across samples) for visual comparability
zscore_rows <- function(mat) {
  t(scale(t(mat)))
}

if_z    <- zscore_rows(if_sub)
serum_z <- zscore_rows(serum_sub)

# Common color breaks (±2 SD)
bk <- seq(-2, 2, length.out = 101)
col_fun <- colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9",
                              "#E0F3F8","#FFFFBF","#FEE090","#FDAE61",
                              "#F46D43","#D73027","#A50026"))(100)

pdf(file.path(out_dir, "heatmap_common_proteins_IF.pdf"), width = 14, height = 5)
pheatmap(if_z,
         color        = col_fun,
         breaks       = bk,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main         = "IF – common proteins (z-score NPX)",
         fontsize_row = 10,
         fontsize_col = 6,
         border_color = NA)
dev.off()

pdf(file.path(out_dir, "heatmap_common_proteins_SERUM.pdf"), width = 14, height = 5)
pheatmap(serum_z,
         color        = col_fun,
         breaks       = bk,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main         = "SERUM – common proteins (z-score NPX)",
         fontsize_row = 10,
         fontsize_col = 6,
         border_color = NA)
dev.off()

# ── 2. Per-protein Pearson correlation (paired samples) ───────────────────────
cor_results <- data.frame(
  protein   = common_proteins,
  r         = NA_real_,
  p_value   = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(common_proteins)) {
  prot <- common_proteins[i]
  ct   <- cor.test(as.numeric(if_sub[prot, ]),
                   as.numeric(serum_sub[prot, ]),
                   method = "pearson")
  cor_results$r[i]       <- ct$estimate
  cor_results$p_value[i] <- ct$p.value
}

cor_results$p_adj <- p.adjust(cor_results$p_value, method = "BH")
cor_results        <- cor_results[order(cor_results$r, decreasing = TRUE), ]

write.table(cor_results,
            file.path(out_dir, "IF_vs_SERUM_common_proteins_correlation.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
print(cor_results)

# ── 3. Correlation bar plot ────────────────────────────────────────────────────
cor_results$protein <- factor(cor_results$protein,
                               levels = cor_results$protein)  # keep sorted order

p_cor <- ggplot(cor_results, aes(x = protein, y = r,
                                  fill = ifelse(p_adj < 0.05, "sig", "ns"))) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(sig = "#D73027", ns = "grey65"),
                    name   = "FDR < 0.05") +
  labs(title = "Pearson r: IF vs SERUM (paired samples)",
       x = NULL, y = "Pearson r") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "barplot_correlation_IF_vs_SERUM.pdf"),
       p_cor, width = 8, height = 4)

# ── 4. Scatter plots per protein ──────────────────────────────────────────────
plot_list <- lapply(common_proteins, function(prot) {
  df <- data.frame(
    IF    = as.numeric(if_sub[prot, ]),
    SERUM = as.numeric(serum_sub[prot, ])
  )
  r_val <- cor_results$r[cor_results$protein == prot]
  p_val <- cor_results$p_adj[cor_results$protein == prot]
  label <- sprintf("r = %.2f\nFDR = %.3f", r_val, p_val)

  ggplot(df, aes(x = IF, y = SERUM)) +
    geom_point(size = 1.5, colour = "#2166AC", alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, colour = "#D73027", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, label = label,
             hjust = -0.1, vjust = 1.3, size = 3) +
    labs(title = prot, x = "IF (NPX)", y = "SERUM (NPX)") +
    theme_bw(base_size = 10)
})

pdf(file.path(out_dir, "scatterplots_IF_vs_SERUM_per_protein.pdf"),
    width = 12, height = 10)
n_col <- 4
n_row <- ceiling(length(plot_list) / n_col)
gridExtra::grid.arrange(grobs = plot_list, ncol = n_col)
dev.off()

# ── 5. Overall IF vs SERUM concordance (all proteins stacked) ─────────────────
if_long    <- melt(as.matrix(if_sub),    varnames = c("protein","sample"),
                   value.name = "NPX_IF")
serum_long <- melt(as.matrix(serum_sub), varnames = c("protein","sample"),
                   value.name = "NPX_SERUM")
joint      <- merge(if_long, serum_long, by = c("protein","sample"))
joint$sample  <- as.integer(as.character(joint$sample))

p_all <- ggplot(joint, aes(x = NPX_IF, y = NPX_SERUM, colour = protein)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(group = protein), linewidth = 0.6) +
  labs(title = "IF vs SERUM expression – all common proteins",
       x = "IF (NPX)", y = "SERUM (NPX)", colour = "Protein") +
  theme_bw(base_size = 11) +
  facet_wrap(~protein, scales = "free", ncol = 4) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "scatter_faceted_IF_vs_SERUM.pdf"),
       p_all, width = 14, height = 12)

cat("\nDone. Output files written to:", out_dir, "\n")
