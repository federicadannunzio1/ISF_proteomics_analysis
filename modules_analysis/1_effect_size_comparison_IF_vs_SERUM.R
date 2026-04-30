rm(list = ls())

# ── Purpose ────────────────────────────────────────────────────────────────────
# Compare the T2D effect size (Cohen's d) for each common protein between
# IF and SERUM. NPX values cannot be compared directly across compartments, but
# Cohen's d (standardized mean difference) is unit-free and comparable.
#
# Positive d  →  protein is UP in T2D
# Negative d  →  protein is DOWN in T2D
#
# Key question: do the same proteins show a stronger T2D signal in IF vs serum?
#
# Outputs:
#   effect_sizes_common_proteins.txt         – table with d, p, FDR per compartment
#   scatter_effect_sizes_IF_vs_SERUM.pdf     – IF d vs SERUM d scatter (quadrant plot)
#   dumbbell_effect_sizes_IF_vs_SERUM.pdf    – dumbbell chart sorted by IF effect
# ──────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(reshape2)
library(ggrepel)

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

cohen_d <- function(group_T2D, group_Control) {
  # Pooled Cohen's d; positive = T2D > Control
  n1 <- length(group_T2D);  n2 <- length(group_Control)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  sp <- sqrt(((n1 - 1) * var(group_T2D) + (n2 - 1) * var(group_Control)) /
               (n1 + n2 - 2))
  if (sp == 0) return(NA_real_)
  (mean(group_T2D) - mean(group_Control)) / sp
}

# ── Load data ──────────────────────────────────────────────────────────────────
if_mat       <- read_mat(file.path(wgcna_dir, "IF/matrix/matrix.txt"))
serum_mat    <- read_mat(file.path(wgcna_dir, "SERUM/matrix/matrix.txt"))
if_traits    <- read_traits(file.path(wgcna_dir, "IF/matrix/Traits.txt"))
serum_traits <- read_traits(file.path(wgcna_dir, "SERUM/matrix/Traits.txt"))

# ── Common proteins ────────────────────────────────────────────────────────────
common_proteins <- c("APOM","ENG","TGFBI","DPP4","MET","DNER",
                     "NCAM1","CHL1","TCN2","PLXNB2","TNXB","FETUB","IL7R")
common_proteins <- common_proteins[common_proteins %in% rownames(if_mat) &
                                     common_proteins %in% rownames(serum_mat)]
cat("Common proteins retained:", length(common_proteins), "\n")

# ── Compute Cohen's d + Wilcoxon p per protein per compartment ─────────────────
compute_effects <- function(mat, traits, proteins, t2d_col = "T2D") {
  shared <- intersect(colnames(mat), rownames(traits))
  mat_s  <- mat[proteins, shared, drop = FALSE]
  t2d    <- traits[shared, t2d_col]

  do.call(rbind, lapply(proteins, function(prot) {
    vals      <- as.numeric(mat_s[prot, ])
    t2d_vals  <- vals[t2d == 1]
    ctrl_vals <- vals[t2d == 0]
    d  <- cohen_d(t2d_vals, ctrl_vals)
    wt <- wilcox.test(t2d_vals, ctrl_vals, exact = FALSE)
    data.frame(protein = prot, cohen_d = d, p_wilcox = wt$p.value,
               n_T2D = length(t2d_vals), n_Control = length(ctrl_vals),
               stringsAsFactors = FALSE)
  }))
}

ef_if    <- compute_effects(if_mat,    if_traits,    common_proteins)
ef_serum <- compute_effects(serum_mat, serum_traits, common_proteins)

ef_if$p_adj    <- p.adjust(ef_if$p_wilcox,    method = "BH")
ef_serum$p_adj <- p.adjust(ef_serum$p_wilcox, method = "BH")

# ── Merge and save ─────────────────────────────────────────────────────────────
ef <- merge(ef_if[, c("protein","cohen_d","p_wilcox","p_adj")],
            ef_serum[, c("protein","cohen_d","p_wilcox","p_adj")],
            by = "protein", suffixes = c("_IF","_SERUM"))

ef$IF_stronger <- abs(ef$cohen_d_IF) > abs(ef$cohen_d_SERUM)
ef <- ef[order(ef$cohen_d_IF), ]

write.table(ef, file.path(out_dir, "effect_sizes_common_proteins.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n── Effect sizes (Cohen's d) ──\n")
print(ef[, c("protein","cohen_d_IF","cohen_d_SERUM","p_adj_IF","p_adj_SERUM","IF_stronger")])
cat("\nProteins with stronger effect in IF:",
    sum(ef$IF_stronger, na.rm = TRUE), "/", nrow(ef), "\n")

# ── Classify significance ──────────────────────────────────────────────────────
ef$sig_class <- with(ef, ifelse(
  p_adj_IF < 0.05 & p_adj_SERUM < 0.05, "Both sig. (FDR<0.05)",
  ifelse(p_adj_IF < 0.05,    "IF only sig.",
  ifelse(p_adj_SERUM < 0.05, "SERUM only sig.", "Neither sig."))))

# ── Plot 1: Scatter plot IF d vs SERUM d (quadrant) ───────────────────────────
lim <- max(abs(c(ef$cohen_d_IF, ef$cohen_d_SERUM)), na.rm = TRUE) * 1.15

p_scatter <- ggplot(ef, aes(x = cohen_d_SERUM, y = cohen_d_IF,
                              colour = sig_class, label = protein)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  # Identity line: points above = stronger effect in IF
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey30", linewidth = 0.5) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 3.2, max.overlaps = 20, show.legend = FALSE,
                  fontface = "bold") +
  scale_colour_manual(
    values = c("Both sig. (FDR<0.05)" = "#D73027",
               "IF only sig."         = "#4575B4",
               "SERUM only sig."      = "#74ADD1",
               "Neither sig."         = "grey60")) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  annotate("text", x = -lim * 0.9, y = lim * 0.9,
           label = "Stronger in IF", size = 3, colour = "grey30", hjust = 0) +
  annotate("text", x = lim * 0.9, y = -lim * 0.9,
           label = "Stronger in SERUM", size = 3, colour = "grey30", hjust = 1) +
  labs(title   = "T2D effect size (Cohen's d): IF vs SERUM",
       subtitle = "Dotted diagonal = equal effect; above diagonal = IF effect is stronger",
       x = "Cohen's d  SERUM  (T2D − Control)",
       y = "Cohen's d  IF  (T2D − Control)",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9))

ggsave(file.path(out_dir, "scatter_effect_sizes_IF_vs_SERUM.pdf"),
       p_scatter, width = 7, height = 7)

# ── Plot 2: Dumbbell chart ─────────────────────────────────────────────────────
ef_long <- melt(ef[, c("protein","cohen_d_IF","cohen_d_SERUM")],
                id.vars = "protein",
                variable.name = "compartment", value.name = "cohen_d")
ef_long$compartment <- gsub("cohen_d_", "", ef_long$compartment)
ef_long$protein     <- factor(ef_long$protein,
                               levels = ef$protein[order(ef$cohen_d_IF)])

p_dumb <- ggplot(ef_long, aes(x = cohen_d, y = protein, colour = compartment)) +
  geom_line(aes(group = protein), colour = "grey70", linewidth = 1) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = c(IF = "#4575B4", SERUM = "#D73027")) +
  labs(title    = "Effect size comparison: IF vs SERUM (common proteins)",
       subtitle  = "T2D vs Control  |  sorted by IF Cohen's d\nBlue = IF  |  Red = SERUM  |  Segment length = discordance between compartments",
       x        = "Cohen's d  (negative = DOWN in T2D)",
       y        = NULL,
       colour   = "Compartment") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "dumbbell_effect_sizes_IF_vs_SERUM.pdf"),
       p_dumb, width = 7, height = 6)

cat("\nDone. Output written to:", out_dir, "\n")
