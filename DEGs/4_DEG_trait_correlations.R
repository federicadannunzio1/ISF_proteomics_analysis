rm(list = ls())

library(corrplot)
library(RColorBrewer)

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
base_path   <- "~/Documents/projects/ISF/ISF_fede"
path_DEG_IF <- file.path(base_path, "code/DEGs/res/IF/DEG.txt")
path_DEG_SE <- file.path(base_path, "code/DEGs/res/SERUM/DEG.txt")

path_mat_IF <- file.path(base_path, "data/matrix/IF_joint.txt")
path_mat_SE <- file.path(base_path, "data/matrix/SERUM_joint.txt")

path_tr_IF  <- file.path(base_path, "data/traits/Traits_IF.txt")
path_tr_SE  <- file.path(base_path, "data/traits/Traits_SERUM.txt")

out_IF      <- file.path(base_path, "code/DEGs/res/IF")
out_SE      <- file.path(base_path, "code/DEGs/res/SERUM")

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────
soglia_pval  <- 0.05
thr_frac_NA  <- 0.90   # exclude traits with >90% NA
thr_frac_zer <- 0.90   # exclude traits with >90% zeros

# ─────────────────────────────────────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────────────────────────────────────
get_corrplot <- function(rho, pval, filepath, soglia_pval) {
  colors <- colorRampPalette(
    c(rev(brewer.pal(9, "Blues")), "white", brewer.pal(9, "Reds"))
  )(200)
  pdf(file = filepath)
  corrplot(rho,
           is.corr      = FALSE,   # must be FALSE for non-square (proteins x traits) matrices
           method       = "color",
           addgrid.col  = "grey80",
           p.mat        = pval,
           sig.level    = soglia_pval,
           insig        = "blank",
           tl.cex       = 0.6,
           cl.cex       = 0.5,
           cl.length    = 5,
           na.label     = "square",
           na.label.col = "white",
           tl.col       = "black",
           col          = colors)
  dev.off()
}

# Replot a corrplot adding extra (non-significant) trait columns as blank.
# Uses res$rho and res$pval (full matrices, guaranteed same dims) — avoids
# any dimension issues from pval_sel.
rebuild_corrplot <- function(res, extra_traits, filepath, soglia_pval) {
  if (length(extra_traits) == 0) return(invisible(NULL))

  all_traits <- c(colnames(res$rho_sel), extra_traits)
  prot_order <- rownames(res$rho_sel)   # alphabetical order from compute_corrplot

  rho_full  <- res$rho[prot_order,  all_traits, drop = FALSE]
  pval_full <- res$pval[prot_order, all_traits, drop = FALSE]

  # Extra trait columns are always blank (p = 1 → non-significant)
  pval_full[, extra_traits] <- 1

  get_corrplot(rho_full, pval_full, filepath, soglia_pval)
}

# Main correlation function (mirroring ISF_giulia compute_corrplot)
compute_corrplot <- function(mat_deg, traits, tag, out_dir, soglia_pval) {

  # Remove traits with too many NA or zeros
  frac_NA  <- apply(traits, 2, function(x) mean(is.na(x)))
  frac_zer <- apply(traits, 2, function(x) mean(x == 0, na.rm = TRUE))
  traits   <- traits[, frac_NA <= thr_frac_NA & frac_zer <= thr_frac_zer,
                     drop = FALSE]

  # Align samples
  common_samples <- intersect(colnames(mat_deg), rownames(traits))
  if (length(common_samples) == 0)
    stop("No common samples between matrix and traits for tag: ", tag)

  data   <- t(mat_deg[, common_samples, drop = FALSE])  # samples × proteins
  traits <- traits[common_samples, , drop = FALSE]

  # Pearson correlation
  # Use cor() directly to guarantee rho, pval, pval_adj all have identical
  # dimensions (proteins × traits). corr.test() returns $p in a non-rectangular
  # format for cross-matrix inputs, causing downstream dimension mismatches.
  rho <- cor(data, traits, use = "pairwise.complete.obs", method = "pearson")

  # Pairwise sample sizes (needed for correct df when NAs are present)
  n_mat <- crossprod(!is.na(data), !is.na(traits))

  # Two-tailed p-value from t-distribution
  t_stat <- rho * sqrt(n_mat - 2) / sqrt(1 - rho^2)
  pval   <- 2 * pt(-abs(t_stat), df = n_mat - 2)
  dimnames(pval) <- dimnames(rho)

  # FDR correction across all protein-trait pairs simultaneously
  pval_adj <- matrix(p.adjust(as.vector(pval), method = "fdr"),
                     nrow = nrow(rho), ncol = ncol(rho),
                     dimnames = dimnames(rho))

  # Select traits with >= 1 significant correlation at soglia_pval
  sig_idx  <- which(pval <= soglia_pval, arr.ind = TRUE)
  if (nrow(sig_idx) == 0) {
    message("No significant correlations found for tag: ", tag)
    sel_cols <- character(0)
  } else {
    sel_cols <- unique(sig_idx[, "col"])
  }

  rho_sel      <- rho[,  sel_cols, drop = FALSE]
  pval_sel     <- pval[, sel_cols, drop = FALSE]
  pval_adj_sel <- pval_adj[, sel_cols, drop = FALSE]

  # Order rows alphabetically (as in ISF_giulia)
  rho_sel      <- rho_sel[order(rownames(rho_sel)), , drop = FALSE]
  pval_sel     <- pval_sel[rownames(rho_sel), , drop = FALSE]
  pval_adj_sel <- pval_adj_sel[rownames(rho_sel), , drop = FALSE]

  # Plot (significant traits only)
  if (ncol(rho_sel) > 0) {
    get_corrplot(rho_sel, pval_sel,
                 file.path(out_dir, paste0("corrplot_", tag, ".pdf")),
                 soglia_pval)
  }

  # Table: significant traits
  if (ncol(rho_sel) > 0) {
    df <- data.frame(cbind(rho_sel, pval_sel, pval_adj_sel))
    colnames(df) <- c(paste0(colnames(rho_sel), "_corr"),
                      paste0(colnames(pval_sel), "_pval"),
                      paste0(colnames(pval_adj_sel), "_pval_adj"))
    df <- df[, order(colnames(df))]
    df <- df[order(rownames(df)), ]
    write.table(df,
                file = file.path(out_dir, paste0("table_", tag, ".txt")),
                quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
  }

  # Table: all traits
  df_all <- data.frame(cbind(rho, pval, pval_adj))
  colnames(df_all) <- c(paste0(colnames(rho), "_corr"),
                        paste0(colnames(pval), "_pval"),
                        paste0(colnames(pval_adj), "_pval_adj"))
  df_all <- df_all[, order(colnames(df_all))]
  df_all <- df_all[order(rownames(df_all)), ]
  write.table(df_all,
              file = file.path(out_dir, paste0("table_all_", tag, ".txt")),
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

  list(rho = rho, pval = pval, pval_adj = pval_adj,
       rho_sel = rho_sel, pval_sel = pval_sel, pval_adj_sel = pval_adj_sel)
}

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────
DEG_IF <- read.delim(path_DEG_IF)$protein
DEG_SE <- read.delim(path_DEG_SE)$protein

mat_IF_full <- read.delim(path_mat_IF, check.names = FALSE, row.names = 1)
mat_SE_full <- read.delim(path_mat_SE, check.names = FALSE, row.names = 1)

# Subset to DEGs (proteins × samples)
# NOTE: Olink NPX data is already on a log2 scale; no further transformation applied.
IF_x <- mat_IF_full[rownames(mat_IF_full) %in% DEG_IF, ]
SE_x <- mat_SE_full[rownames(mat_SE_full) %in% DEG_SE, ]

# Load traits (samples × traits, row.names = sample IDs)
traits_IF <- read.delim(path_tr_IF,  check.names = FALSE, row.names = 1,
                        na.strings = c("NA", " ", "", "n/a", "N/A"))
traits_SE <- read.delim(path_tr_SE,  check.names = FALSE, row.names = 1,
                        na.strings = c("NA", " ", "", "n/a", "N/A"))

# Remove T2D column (used for stratification, not as a continuous trait)
traits_IF <- traits_IF[, colnames(traits_IF) != "T2D", drop = FALSE]
traits_SE <- traits_SE[, colnames(traits_SE) != "T2D", drop = FALSE]

# Extract T2D status from matrix column alignment using original traits file
t2d_IF <- read.delim(path_tr_IF, check.names = FALSE, row.names = 1,
                     na.strings = c("NA", " ", "", "n/a", "N/A"))[, "T2D"]
names(t2d_IF) <- rownames(read.delim(path_tr_IF, check.names = FALSE,
                                     row.names = 1))

t2d_SE <- read.delim(path_tr_SE, check.names = FALSE, row.names = 1,
                     na.strings = c("NA", " ", "", "n/a", "N/A"))[, "T2D"]
names(t2d_SE) <- rownames(read.delim(path_tr_SE, check.names = FALSE,
                                     row.names = 1))

# Split samples by T2D status
IF_cases    <- names(t2d_IF)[t2d_IF == 1]
IF_controls <- names(t2d_IF)[t2d_IF == 0]
SE_cases    <- names(t2d_SE)[t2d_SE == 1]
SE_controls <- names(t2d_SE)[t2d_SE == 0]

# ─────────────────────────────────────────────────────────────────────────────
# IF correlations
# ─────────────────────────────────────────────────────────────────────────────
message("--- IF: T2D cases ---")
res_IF_T2D <- compute_corrplot(
  mat_deg    = IF_x[, colnames(IF_x) %in% IF_cases,    drop = FALSE],
  traits     = traits_IF[rownames(traits_IF) %in% IF_cases, , drop = FALSE],
  tag        = "T2D",
  out_dir    = out_IF,
  soglia_pval = soglia_pval
)

message("--- IF: controls ---")
res_IF_ctrl <- compute_corrplot(
  mat_deg    = IF_x[, colnames(IF_x) %in% IF_controls, drop = FALSE],
  traits     = traits_IF[rownames(traits_IF) %in% IF_controls, , drop = FALSE],
  tag        = "control",
  out_dir    = out_IF,
  soglia_pval = soglia_pval
)

# Rebuild plots adding columns unique to the other group (as in ISF_giulia)
extra_IF_T2D  <- setdiff(colnames(res_IF_ctrl$rho_sel), colnames(res_IF_T2D$rho_sel))
extra_IF_ctrl <- setdiff(colnames(res_IF_T2D$rho_sel),  colnames(res_IF_ctrl$rho_sel))

rebuild_corrplot(res_IF_T2D,  extra_IF_T2D,
                 file.path(out_IF, "corrplot_T2D_extended.pdf"),    soglia_pval)
rebuild_corrplot(res_IF_ctrl, extra_IF_ctrl,
                 file.path(out_IF, "corrplot_control_extended.pdf"), soglia_pval)

# ─────────────────────────────────────────────────────────────────────────────
# SERUM correlations
# ─────────────────────────────────────────────────────────────────────────────
message("--- SERUM: T2D cases ---")
res_SE_T2D <- compute_corrplot(
  mat_deg    = SE_x[, colnames(SE_x) %in% SE_cases,    drop = FALSE],
  traits     = traits_SE[rownames(traits_SE) %in% SE_cases, , drop = FALSE],
  tag        = "T2D",
  out_dir    = out_SE,
  soglia_pval = soglia_pval
)

message("--- SERUM: controls ---")
res_SE_ctrl <- compute_corrplot(
  mat_deg    = SE_x[, colnames(SE_x) %in% SE_controls, drop = FALSE],
  traits     = traits_SE[rownames(traits_SE) %in% SE_controls, , drop = FALSE],
  tag        = "control",
  out_dir    = out_SE,
  soglia_pval = soglia_pval
)

# Rebuild plots
extra_SE_T2D  <- setdiff(colnames(res_SE_ctrl$rho_sel), colnames(res_SE_T2D$rho_sel))
extra_SE_ctrl <- setdiff(colnames(res_SE_T2D$rho_sel),  colnames(res_SE_ctrl$rho_sel))

rebuild_corrplot(res_SE_T2D,  extra_SE_T2D,
                 file.path(out_SE, "corrplot_T2D_extended.pdf"),    soglia_pval)
rebuild_corrplot(res_SE_ctrl, extra_SE_ctrl,
                 file.path(out_SE, "corrplot_control_extended.pdf"), soglia_pval)

message("Done.")

# ─────────────────────────────────────────────────────────────────────────────
# Combined figure: IF (top) | SERUM (bottom), T2D | Control side by side
# Same clinical-trait order on x-axis; proteins labelled on y-axis
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2)
library(patchwork)

# Unified x-axis: traits significant in at least one group across IF + SERUM
# (same logic as the "extended" corrplots, but spanning all 4 groups)
sig_traits_all <- sort(unique(c(
  colnames(res_IF_T2D$rho_sel),  colnames(res_IF_ctrl$rho_sel),
  colnames(res_SE_T2D$rho_sel),  colnames(res_SE_ctrl$rho_sel)
)))

if (length(sig_traits_all) == 0) {
  warning("No significant traits in any group — falling back to all traits for combined figure.")
  sig_traits_all <- sort(unique(c(
    colnames(res_IF_T2D$rho), colnames(res_IF_ctrl$rho),
    colnames(res_SE_T2D$rho), colnames(res_SE_ctrl$rho)
  )))
}

# Protein orders per fluid (alphabetical; union in case T2D/ctrl differ due to NA filter)
prot_order_IF <- sort(unique(c(rownames(res_IF_T2D$rho), rownames(res_IF_ctrl$rho))))
prot_order_SE <- sort(unique(c(rownames(res_SE_T2D$rho), rownames(res_SE_ctrl$rho))))

# Color palette matching get_corrplot (reversed Blues -> white -> Reds)
pal_colors <- colorRampPalette(
  c(rev(brewer.pal(9, "Blues")), "white", brewer.pal(9, "Reds"))
)(200)
low_col  <- pal_colors[1]
high_col <- pal_colors[200]

# Helper: build one ggplot tile panel
# rho, pval: full matrices from res$rho / res$pval (proteins x traits)
make_gg_tile <- function(rho, pval, title, soglia_pval, trait_order, prot_order) {

  # Expand to unified trait_order and prot_order (fill absent entries with NA / p=1)
  rho_exp  <- matrix(NA_real_, nrow = length(prot_order), ncol = length(trait_order),
                     dimnames = list(prot_order, trait_order))
  pval_exp <- matrix(1,        nrow = length(prot_order), ncol = length(trait_order),
                     dimnames = list(prot_order, trait_order))

  shared_prot  <- intersect(prot_order,  rownames(rho))
  shared_trait <- intersect(trait_order, colnames(rho))
  rho_exp[ shared_prot, shared_trait] <- rho[ shared_prot, shared_trait]
  pval_exp[shared_prot, shared_trait] <- pval[shared_prot, shared_trait]

  # Mask non-significant cells (show as NA = white tile)
  rho_plot <- rho_exp
  rho_plot[pval_exp > soglia_pval] <- NA

  # Long format using base R
  df <- data.frame(
    protein = rep(prot_order, times = length(trait_order)),
    trait   = rep(trait_order, each  = length(prot_order)),
    rho     = as.vector(rho_plot),
    stringsAsFactors = FALSE
  )
  # Reverse protein factor so alphabetical order reads top-to-bottom on y-axis
  df$protein <- factor(df$protein, levels = rev(prot_order))
  df$trait   <- factor(df$trait,   levels = trait_order)

  ggplot(df, aes(x = trait, y = protein, fill = rho)) +
    geom_tile(color = "grey80", linewidth = 0.3) +
    scale_fill_gradientn(
      colors   = pal_colors,
      limits   = c(-1, 1),
      na.value = "white",
      name     = "r",
      guide    = guide_colorbar(barheight = 8, barwidth = 0.8)
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y  = element_text(size = 7),
      plot.title   = element_text(size = 10, face = "bold"),
      panel.grid   = element_blank()
    )
}

# Build 4 panels
p_IF_T2D  <- make_gg_tile(res_IF_T2D$rho,  res_IF_T2D$pval,  "IF - T2D",
                           soglia_pval, sig_traits_all, prot_order_IF)
p_IF_ctrl <- make_gg_tile(res_IF_ctrl$rho,  res_IF_ctrl$pval, "IF - Control",
                           soglia_pval, sig_traits_all, prot_order_IF)
p_SE_T2D  <- make_gg_tile(res_SE_T2D$rho,  res_SE_T2D$pval,  "SERUM - T2D",
                           soglia_pval, sig_traits_all, prot_order_SE)
p_SE_ctrl <- make_gg_tile(res_SE_ctrl$rho,  res_SE_ctrl$pval, "SERUM - Control",
                           soglia_pval, sig_traits_all, prot_order_SE)

# Assemble: top row = IF, bottom row = SERUM; legend shared
fig_combined <-
  (p_IF_T2D | p_IF_ctrl) /
  (p_SE_T2D | p_SE_ctrl) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

out_combined <- file.path(base_path, "code/DEGs/res")
ggsave(file.path(out_combined, "combined_corrplot.pdf"),
       plot = fig_combined, width = 14, height = 10)
ggsave(file.path(out_combined, "combined_corrplot.png"),
       plot = fig_combined, width = 14, height = 10, dpi = 300)

message("Combined figure saved in: ", out_combined)
