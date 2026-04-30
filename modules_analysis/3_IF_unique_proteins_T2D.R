rm(list = ls())

# ── Purpose ────────────────────────────────────────────────────────────────────
# Analyse the 36 proteins that are EXCLUSIVE to the IF turquoise WGCNA module
# (i.e. present in the IF T2D-associated module but NOT in the SERUM brown module).
#
# These proteins represent the IF-specific proteome signal for T2D and are
# the most biologically interesting from the perspective of "what IF captures
# that serum does not".
#
# Tests performed:
#   - Wilcoxon T2D vs Control for each protein in IF
#   - FDR correction (Benjamini-Hochberg)
#
# Outputs:
#   IF_unique_proteins_T2D_results.txt        – full results table
#   volcano_IF_unique_proteins.pdf            – volcano plot
#   boxplots_IF_unique_proteins_sig.pdf       – boxplots for FDR < 0.05 proteins
# ──────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)

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
if_mat    <- read_mat(file.path(wgcna_dir, "IF/matrix/matrix.txt"))
if_traits <- read_traits(file.path(wgcna_dir, "IF/matrix/Traits.txt"))

# ── IF-unique module proteins (from venn/summary_table) ───────────────────────
if_unique_all <- c("GP1BA","SAA4","NOTCH1","THBS4","SELL","IGFBP3","CRTAC1",
                   "ANGPTL3","TIE1","F11","FCN2","MEGF9","CNDP1","SPARCL1",
                   "C1QTNF1","CFHR5","SERPINA5","SERPINA7","MBL2","ST6GAL1",
                   "PCOLCE","CR2","AOC3","FCGR3B","PROC","NID1","VCAM1",
                   "ICAM3","F7","IL-17A","CA4","TRAIL","ICAM1","ANG","C2","CES1")

if_unique <- if_unique_all[if_unique_all %in% rownames(if_mat)]
cat("IF-unique proteins listed:", length(if_unique_all), "\n")
cat("Found in matrix:          ", length(if_unique), "\n")
if (length(setdiff(if_unique_all, if_unique)) > 0)
  cat("Missing from matrix:", paste(setdiff(if_unique_all, if_unique), collapse=", "), "\n")

# ── Align with traits ─────────────────────────────────────────────────────────
shared    <- intersect(colnames(if_mat), rownames(if_traits))
mat_sub   <- if_mat[if_unique, shared, drop = FALSE]
t2d_vec   <- if_traits[shared, "T2D"]
t2d_label <- factor(ifelse(t2d_vec == 1, "T2D", "Control"), levels = c("Control","T2D"))

cat("IF samples used:", length(shared),
    " | T2D:", sum(t2d_vec == 1), " | Control:", sum(t2d_vec == 0), "\n")

# ── Wilcoxon + delta median per protein ──────────────────────────────────────
results <- do.call(rbind, lapply(if_unique, function(prot) {
  vals      <- as.numeric(mat_sub[prot, ])
  t2d_vals  <- vals[t2d_label == "T2D"]
  ctrl_vals <- vals[t2d_label == "Control"]
  wt        <- wilcox.test(t2d_vals, ctrl_vals, exact = FALSE)
  delta     <- median(t2d_vals, na.rm = TRUE) - median(ctrl_vals, na.rm = TRUE)
  data.frame(
    protein      = prot,
    delta_median = round(delta, 4),
    direction    = ifelse(delta > 0, "UP in T2D", "DOWN in T2D"),
    p_wilcox     = wt$p.value,
    stringsAsFactors = FALSE)
}))

results$p_adj <- p.adjust(results$p_wilcox, method = "BH")
results <- results[order(results$p_adj), ]

write.table(results, file.path(out_dir, "IF_unique_proteins_T2D_results.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n── IF-unique proteins: T2D vs Control ──\n")
print(results)
sig_prots <- results$protein[results$p_adj < 0.05]
cat("\nSignificant at FDR < 0.05:", length(sig_prots),
    if (length(sig_prots) > 0) paste("–", paste(sig_prots, collapse=", ")) else "", "\n")

nom_prots <- results$protein[results$p_wilcox < 0.05 & results$p_adj >= 0.05]
cat("Nominally significant (p<0.05, FDR≥0.05):", length(nom_prots),
    if (length(nom_prots) > 0) paste("–", paste(nom_prots, collapse=", ")) else "", "\n")

# ── Volcano plot ───────────────────────────────────────────────────────────────
# FDR threshold line at the highest p_wilcox that still has p_adj < 0.05
fdr_cutoff_p <- if (any(results$p_adj < 0.05))
  max(results$p_wilcox[results$p_adj < 0.05], na.rm = TRUE) else NA

results$neglog10p <- -log10(results$p_wilcox)
results$label_col <- with(results, ifelse(
  p_adj < 0.05 & delta_median < 0, "Down in T2D (FDR<0.05)",
  ifelse(p_adj < 0.05 & delta_median > 0, "Up in T2D (FDR<0.05)",
  ifelse(p_wilcox < 0.05, "Nominally sig. (p<0.05)", "ns"))))

results$plot_label <- with(results,
  ifelse(p_adj < 0.05 | p_wilcox < 0.01, protein, ""))

p_volc <- ggplot(results, aes(x = delta_median, y = neglog10p,
                               colour = label_col, label = plot_label)) +
  {if (!is.na(fdr_cutoff_p))
    geom_hline(yintercept = -log10(fdr_cutoff_p),
               linetype = "dashed", colour = "grey40", linewidth = 0.5)} +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 25, show.legend = FALSE, fontface = "bold") +
  scale_colour_manual(
    values = c("Down in T2D (FDR<0.05)"    = "#4575B4",
               "Up in T2D (FDR<0.05)"      = "#D73027",
               "Nominally sig. (p<0.05)"   = "#FDAE61",
               "ns"                        = "grey60")) +
  labs(
    title    = "IF-unique module proteins: T2D vs Control (IF compartment)",
    subtitle = "Proteins in IF turquoise WGCNA module only (absent from SERUM brown module)\nDashed line = FDR 5% threshold",
    x        = "Δ median NPX  (T2D − Control)",
    y        = "−log₁₀(p-value)",
    colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "volcano_IF_unique_proteins.pdf"),
       p_volc, width = 8, height = 7)

# ── Boxplots for significant or top-10 proteins ──────────────────────────────
# Show FDR<0.05 proteins; if none, show top 10 by p_wilcox
plot_prots <- if (length(sig_prots) > 0) sig_prots else head(results$protein, 10)
cat("\nGenerating boxplots for:", paste(plot_prots, collapse=", "), "\n")

all_long <- melt(as.matrix(mat_sub[plot_prots, , drop = FALSE]),
                 varnames = c("protein","sample"), value.name = "NPX")
all_long$T2D_label <- t2d_label[as.character(all_long$sample)]
all_long$protein   <- factor(all_long$protein, levels = plot_prots)

ann_sig <- results[results$protein %in% plot_prots, ]
ann_df  <- data.frame(
  protein = ann_sig$protein,
  label   = sapply(ann_sig$p_adj, function(p)
              if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*"
              else sprintf("p=%.3f", ann_sig$p_wilcox[ann_sig$p_adj == p][1])),
  y_pos   = apply(mat_sub[plot_prots, , drop = FALSE], 1,
                  function(x) max(x, na.rm = TRUE) * 1.05)
)
ann_df$protein <- factor(ann_df$protein, levels = plot_prots)

n_col <- min(4, length(plot_prots))
n_row <- ceiling(length(plot_prots) / n_col)

p_box <- ggplot(all_long, aes(x = T2D_label, y = NPX, fill = T2D_label)) +
  geom_boxplot(outlier.size = 0.8, width = 0.55, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  geom_text(data = ann_df, aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4.5) +
  scale_fill_manual(values = c(Control = "#6BAED6", T2D = "#FB6A4A")) +
  facet_wrap(~protein, scales = "free_y", ncol = n_col) +
  labs(
    title    = "IF-unique module proteins: T2D vs Control",
    subtitle = if (length(sig_prots) > 0) "FDR < 0.05" else "Top 10 by p-value (no FDR<0.05)",
    x = NULL, y = "NPX") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold"))

ggsave(file.path(out_dir, "boxplots_IF_unique_proteins_sig.pdf"),
       p_box, width = 4 * n_col, height = 4 * n_row)

cat("\nDone. Output written to:", out_dir, "\n")
