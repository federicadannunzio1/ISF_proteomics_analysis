rm(list = ls())

# ── Purpose ────────────────────────────────────────────────────────────────────
# For each paired subject (same individual has both IF and SERUM), compute the
# within-subject z-score difference:
#
#   Δ = z_IF  −  z_SERUM   (z-scored within each compartment separately)
#
# Δ > 0  →  the protein is relatively HIGHER in IF than in SERUM for that subject
# Δ < 0  →  the protein is relatively LOWER  in IF than in SERUM
#
# Then test whether Δ differs between T2D and Control using Wilcoxon.
# A significant result means T2D specifically alters the IF/SERUM balance
# of that protein — i.e., the protein is NOT simply lower everywhere in T2D,
# but is differentially affected in IF relative to serum.
#
# Why z-score first?
#   NPX values are relative within each Olink panel and CANNOT be compared
#   directly across IF and SERUM. Z-scoring within each compartment (across
#   all subjects) puts both compartments on a mean-0, SD-1 scale, making
#   the subtraction meaningful.
#
# Outputs:
#   paired_zdiff_T2D_vs_control.txt        – Wilcoxon test results per protein
#   boxplot_zdiff_per_protein.pdf          – Δ boxplots per protein, T2D vs Control
#   heatmap_zdiff_IF_minus_SERUM.pdf       – heatmap of Δ matrix, sorted by T2D status
# ──────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(reshape2)
library(pheatmap)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir  <- "~/Documents/projects/ISF/ISF_fede"
wgcna_dir <- file.path(base_dir, "code/WGCNA_paola/code/project/ISF/dataset")
out_dir   <- file.path(base_dir, "code/modules_analysis")

# ── Helpers ────────────────────────────────────────────────────────────────────
read_mat <- function(path) {
  m <- read.table(path, sep = "\t", header = TRUE, quote = "",
                  check.names = FALSE, row.names = 1)
  rownames(m) <- gsub('"', '', rownames(m))
  colnames(m) <- gsub('"', '', colnames(m))
  m
}

read_traits <- function(path) {
  t <- read.table(path, sep = "\t", header = TRUE, quote = "",
                  check.names = FALSE, row.names = 1)
  rownames(t) <- gsub('"', '', rownames(t))
  t
}

# ── Load data ──────────────────────────────────────────────────────────────────
if_mat       <- read_mat(file.path(wgcna_dir, "IF/matrix/matrix.txt"))
serum_mat    <- read_mat(file.path(wgcna_dir, "SERUM/matrix/matrix.txt"))
if_traits    <- read_traits(file.path(wgcna_dir, "IF/matrix/Traits.txt"))

# ── Common proteins ────────────────────────────────────────────────────────────
common_proteins <- c("APOM","ENG","TGFBI","DPP4","MET","DNER",
                     "NCAM1","CHL1","TCN2","PLXNB2","TNXB","FETUB","IL7R")
common_proteins <- common_proteins[common_proteins %in% rownames(if_mat) &
                                     common_proteins %in% rownames(serum_mat)]

# ── Align PAIRED samples ───────────────────────────────────────────────────────
# Each subject has one IF sample ("IF X") and one SERUM sample ("Serum X")
# with the same numeric ID X
if_ids    <- as.integer(gsub("^IF\\s+",    "", colnames(if_mat)))
serum_ids <- as.integer(gsub("^Serum\\s+", "", colnames(serum_mat)))
shared_ids <- sort(intersect(if_ids, serum_ids))
cat("Paired subjects:", length(shared_ids), "\n")

if_sub    <- if_mat[common_proteins,    match(shared_ids, if_ids),    drop = FALSE]
serum_sub <- serum_mat[common_proteins, match(shared_ids, serum_ids), drop = FALSE]
colnames(if_sub)    <- as.character(shared_ids)
colnames(serum_sub) <- as.character(shared_ids)

# ── Z-score within each compartment (across paired subjects only) ──────────────
# Result: for each protein, IF and SERUM are both on a N(0,1) scale
zscore_rows <- function(mat) t(scale(t(mat)))
if_z    <- zscore_rows(if_sub)    # proteins × subjects, z-scored within IF
serum_z <- zscore_rows(serum_sub) # proteins × subjects, z-scored within SERUM

# ── Within-subject difference: Δ = z_IF − z_SERUM ────────────────────────────
delta_z <- if_z - serum_z   # proteins × subjects

# ── T2D status for paired subjects ────────────────────────────────────────────
if_sample_names <- colnames(if_mat)[match(shared_ids, if_ids)]
t2d_vec   <- setNames(if_traits[if_sample_names, "T2D"], as.character(shared_ids))
t2d_label <- factor(ifelse(t2d_vec == 1, "T2D", "Control"), levels = c("Control","T2D"))

# ── Wilcoxon test: does Δ differ between T2D and Control? ─────────────────────
zdiff_res <- do.call(rbind, lapply(common_proteins, function(prot) {
  vals      <- delta_z[prot, ]
  t2d_vals  <- vals[t2d_label == "T2D"]
  ctrl_vals <- vals[t2d_label == "Control"]
  wt <- wilcox.test(t2d_vals, ctrl_vals, exact = FALSE)
  data.frame(
    protein            = prot,
    mean_delta_T2D     = round(mean(t2d_vals,  na.rm = TRUE), 4),
    mean_delta_Control = round(mean(ctrl_vals, na.rm = TRUE), 4),
    # A negative diff means T2D subjects have lower z_IF relative to z_SERUM
    # (i.e. the protein drops MORE in IF than in SERUM in T2D)
    diff_T2D_minus_Control = round(mean(t2d_vals, na.rm=TRUE) - mean(ctrl_vals, na.rm=TRUE), 4),
    p_wilcox = wt$p.value,
    stringsAsFactors = FALSE)
}))

zdiff_res$p_adj <- p.adjust(zdiff_res$p_wilcox, method = "BH")
zdiff_res <- zdiff_res[order(zdiff_res$p_adj), ]

write.table(zdiff_res, file.path(out_dir, "paired_zdiff_T2D_vs_control.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n── Paired z-diff results (Δ = z_IF − z_SERUM; T2D vs Control) ──\n")
print(zdiff_res)

# ── Plot 1: Boxplots of Δ per protein ─────────────────────────────────────────
delta_long            <- melt(t(delta_z), varnames = c("sample","protein"), value.name = "delta_z")
delta_long$sample     <- as.character(delta_long$sample)
delta_long$T2D_label  <- t2d_label[delta_long$sample]
delta_long$protein    <- factor(delta_long$protein, levels = zdiff_res$protein)

# Significance labels
sig_ann <- data.frame(
  protein = zdiff_res$protein,
  label   = sapply(zdiff_res$p_adj, function(p)
              if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"),
  y_pos   = apply(delta_z, 1, function(x) max(x, na.rm = TRUE) + 0.15)
)
sig_ann$protein <- factor(sig_ann$protein, levels = zdiff_res$protein)

p_box <- ggplot(delta_long, aes(x = T2D_label, y = delta_z, fill = T2D_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.8, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.45) +
  geom_text(data = sig_ann, aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4.5) +
  scale_fill_manual(values = c(Control = "#6BAED6", T2D = "#FB6A4A")) +
  facet_wrap(~protein, scales = "free_y", ncol = 4) +
  labs(
    title    = "Within-subject z-score difference: IF − SERUM",
    subtitle = paste0("Δ = z_IF − z_SERUM per subject (z-scored within each compartment)\n",
                      "Δ < 0 in T2D → protein is relatively MORE depleted in IF than in SERUM in T2D\n",
                      "Ordered by FDR"),
    x = NULL, y = "Δ z-score (z_IF − z_SERUM)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold"),
        plot.subtitle   = element_text(size = 8, colour = "grey30"))

ggsave(file.path(out_dir, "boxplot_zdiff_per_protein.pdf"),
       p_box, width = 14, height = 14)

# ── Plot 2: Heatmap of Δ matrix sorted by T2D status ─────────────────────────
sample_order   <- order(t2d_label)
delta_ordered  <- delta_z[, sample_order, drop = FALSE]
ann_col        <- data.frame(Group = t2d_label[sample_order])
rownames(ann_col) <- colnames(delta_ordered)

# Order rows by mean difference (T2D - Control)
row_order <- order(zdiff_res$diff_T2D_minus_Control)
delta_ordered <- delta_ordered[zdiff_res$protein[row_order], , drop = FALSE]

pdf(file.path(out_dir, "heatmap_zdiff_IF_minus_SERUM.pdf"), width = 14, height = 5)
pheatmap(delta_ordered,
         annotation_col    = ann_col,
         annotation_colors = list(Group = c(Control = "#6BAED6", T2D = "#FB6A4A")),
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         main              = "z_IF − z_SERUM per subject  (common proteins)",
         fontsize_row      = 10,
         fontsize_col      = 6,
         border_color      = NA,
         color = colorRampPalette(c("#313695","#FFFFBF","#A50026"))(100))
dev.off()

cat("\nDone. Output written to:", out_dir, "\n")
