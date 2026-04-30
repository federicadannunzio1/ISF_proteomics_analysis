rm(list = ls())

library(ggplot2)
library(reshape2)
library(gridExtra)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir  <- "~/Documents/projects/ISF/ISF_fede"
wgcna_dir <- file.path(base_dir, "code/WGCNA_paola/code/project/ISF/dataset")
out_dir   <- file.path(base_dir, "code/modules_analysis")

# ── Load matrices ──────────────────────────────────────────────────────────────
read_mat <- function(path) {
  m <- read.table(path, sep = "\t", header = TRUE, quote = "",
                  check.names = FALSE, row.names = 1)
  rownames(m) <- gsub('"', '', rownames(m))
  colnames(m) <- gsub('"', '', colnames(m))
  m
}

if_mat    <- read_mat(file.path(wgcna_dir, "IF/matrix/matrix.txt"))
serum_mat <- read_mat(file.path(wgcna_dir, "SERUM/matrix/matrix.txt"))

# ── Load Traits → T2D column ───────────────────────────────────────────────────
read_traits <- function(path) {
  t <- read.table(path, sep = "\t", header = TRUE, quote = "",
                  check.names = FALSE, row.names = 1)
  rownames(t) <- gsub('"', '', rownames(t))
  t
}

if_traits    <- read_traits(file.path(wgcna_dir, "IF/matrix/Traits.txt"))
serum_traits <- read_traits(file.path(wgcna_dir, "SERUM/matrix/Traits.txt"))

# ── Common proteins ────────────────────────────────────────────────────────────
common_proteins <- c("APOM","ENG","TGFBI","DPP4","MET","DNER",
                     "NCAM1","CHL1","TCN2","PLXNB2","TNXB","FETUB","IL7R")
common_proteins <- common_proteins[common_proteins %in% rownames(if_mat) &
                                     common_proteins %in% rownames(serum_mat)]

# ── Build long data frame: compartment + T2D status ───────────────────────────
make_long <- function(mat, traits, compartment, id_prefix) {
  shared_samples <- intersect(colnames(mat), rownames(traits))
  mat_sub <- mat[common_proteins, shared_samples, drop = FALSE]

  # Named vector with explicit names preserved
  t2d_vec <- setNames(traits[shared_samples, "T2D"], shared_samples)

  df <- melt(as.matrix(mat_sub), varnames = c("protein", "sample"), value.name = "NPX")
  df$sample      <- as.character(df$sample)
  df$T2D         <- t2d_vec[df$sample]
  df$compartment <- compartment
  df$sample_id   <- as.integer(gsub(paste0('^', id_prefix, '\\s+'), '', df$sample))
  df
}

if_long    <- make_long(if_mat,    if_traits,    "IF",    "IF")
serum_long <- make_long(serum_mat, serum_traits, "SERUM", "Serum")

all_long <- rbind(if_long, serum_long)
all_long$T2D_label    <- ifelse(all_long$T2D == 1, "T2D", "Control")
all_long$T2D_label    <- factor(all_long$T2D_label, levels = c("Control", "T2D"))
all_long$compartment  <- factor(all_long$compartment, levels = c("IF", "SERUM"))

# ── Compute direction (median T2D − median Control) per protein per compartment
dir_table <- do.call(rbind, lapply(common_proteins, function(prot) {
  sub <- all_long[all_long$protein == prot, ]
  do.call(rbind, lapply(c("IF","SERUM"), function(comp) {
    s      <- sub[sub$compartment == comp, ]
    med_c  <- median(s$NPX[s$T2D_label == "Control"], na.rm = TRUE)
    med_t  <- median(s$NPX[s$T2D_label == "T2D"],     na.rm = TRUE)
    wt     <- wilcox.test(NPX ~ T2D_label, data = s, exact = FALSE)
    data.frame(protein     = prot,
               compartment = comp,
               delta_median = med_t - med_c,
               direction    = ifelse(med_t > med_c, "UP in T2D", "DOWN in T2D"),
               p_wilcox     = wt$p.value,
               stringsAsFactors = FALSE)
  }))
}))

# FDR across all tests
dir_table$p_adj <- p.adjust(dir_table$p_wilcox, method = "BH")

# ── Flag discordant proteins ────────────────────────────────────────────────────
wide_dir <- reshape(dir_table[, c("protein","compartment","direction")],
                    idvar = "protein", timevar = "compartment", direction = "wide")
colnames(wide_dir) <- gsub("direction\\.", "", colnames(wide_dir))
wide_dir$discordant <- wide_dir$IF != wide_dir$SERUM

discordant_proteins <- wide_dir$protein[wide_dir$discordant]
concordant_proteins <- wide_dir$protein[!wide_dir$discordant]

cat("Discordant direction (IF vs SERUM) proteins:\n")
print(discordant_proteins)
cat("\nConcordant direction proteins:\n")
print(concordant_proteins)

# Save direction summary
write.table(merge(wide_dir,
                  reshape(dir_table[, c("protein","compartment","delta_median","p_adj")],
                          idvar="protein", timevar="compartment", direction="wide"),
                  by="protein"),
            file.path(out_dir, "direction_summary_T2D_vs_control.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── Helper: single boxplot for one protein ────────────────────────────────────
make_boxplot <- function(prot, highlight_discordant = FALSE) {
  sub    <- all_long[all_long$protein == prot, ]
  d_info <- dir_table[dir_table$protein == prot, ]

  # Build significance labels per facet
  sig_labels <- sapply(c("IF","SERUM"), function(comp) {
    p <- d_info$p_adj[d_info$compartment == comp]
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  })
  ann_df <- data.frame(compartment = factor(c("IF","SERUM"), levels = c("IF","SERUM")),
                       label       = sig_labels,
                       y_pos       = max(sub$NPX, na.rm = TRUE) * 1.05)

  fill_vals  <- c(Control = "#6BAED6", T2D = "#FB6A4A")
  title_col  <- if (highlight_discordant && prot %in% discordant_proteins) "#C0392B" else "black"

  ggplot(sub, aes(x = T2D_label, y = NPX, fill = T2D_label)) +
    geom_boxplot(outlier.size = 0.8, width = 0.55) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    scale_fill_manual(values = fill_vals) +
    geom_text(data = ann_df, aes(x = 1.5, y = y_pos, label = label),
              inherit.aes = FALSE, size = 4) +
    facet_wrap(~compartment, nrow = 1) +
    labs(title = prot, x = NULL, y = "NPX") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          strip.text      = element_text(face = "bold"),
          plot.title      = element_text(colour = title_col, face = "bold"))
}

# ── Plot 1: DISCORDANT proteins only ──────────────────────────────────────────
if (length(discordant_proteins) > 0) {
  plots_disc <- lapply(discordant_proteins, make_boxplot, highlight_discordant = TRUE)
  n_col <- min(3, length(plots_disc))
  n_row <- ceiling(length(plots_disc) / n_col)

  pdf(file.path(out_dir, "boxplots_DISCORDANT_proteins_T2D_vs_control.pdf"),
      width = 4 * n_col, height = 4 * n_row)
  grid.arrange(grobs = plots_disc, ncol = n_col,
               top = "Proteins with DISCORDANT direction in IF vs SERUM (T2D vs Control)")
  dev.off()
  cat("Discordant boxplots saved.\n")
} else {
  cat("No discordant proteins found.\n")
}

# ── Plot 2: ALL 13 proteins together for reference ────────────────────────────
all_plots <- lapply(common_proteins, make_boxplot, highlight_discordant = TRUE)
pdf(file.path(out_dir, "boxplots_ALL_common_proteins_T2D_vs_control.pdf"),
    width = 14, height = 18)
grid.arrange(grobs = all_plots, ncol = 3,
             top = "All common proteins: T2D vs Control in IF and SERUM\n(red title = discordant direction)")
dev.off()

cat("\nDone. Output files written to:", out_dir, "\n")
print(dir_table[order(dir_table$protein, dir_table$compartment), ])
