rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# ==========================================
# PATHS
# ==========================================
base    <- "~/Documents/projects/ISF/ISF_fede"
dirOut  <- file.path(base, "code", "supplementary_analysis", "output")
deg_dir <- file.path(base, "code", "DEGs", "res")

mat_IF <- read.delim(file.path(base, "data/matrix/IF_joint.txt"),
                     check.names = FALSE, row.names = 1)
mat_SE <- read.delim(file.path(base, "data/matrix/SERUM_joint.txt"),
                     check.names = FALSE, row.names = 1)

tr_IF <- read.delim(file.path(base, "data/traits/Traits_IF.txt"),
                    check.names = FALSE, row.names = 1,
                    na.strings = c("NA", " ", "", "n/a", "N/A"))
tr_SE <- read.delim(file.path(base, "data/traits/Traits_SERUM.txt"),
                    check.names = FALSE, row.names = 1,
                    na.strings = c("NA", " ", "", "n/a", "N/A"))

# ==========================================
# DEG LISTS
# ==========================================
deg_IF_all <- read.delim(file.path(deg_dir, "IF",    "DEG_res", "DEG.txt"))$protein
deg_SE_all <- read.delim(file.path(deg_dir, "SERUM", "DEG_res", "DEG.txt"))$protein

if_unique    <- deg_IF_all[!deg_IF_all %in% deg_SE_all]
serum_unique <- deg_SE_all[!deg_SE_all %in% deg_IF_all]
shared       <- intersect(deg_IF_all, deg_SE_all)

# ==========================================
# PARAMETERS
# ==========================================
pval_thr  <- 0.05
thr_na    <- 0.90
thr_zeros <- 0.90

# Remove T2D column from traits (used for stratification only)
tr_IF <- tr_IF[, colnames(tr_IF) != "T2D", drop = FALSE]
tr_SE <- tr_SE[, colnames(tr_SE) != "T2D", drop = FALSE]

# Remove traits with too many NA or zeros
clean_traits <- function(tr) {
  fna  <- apply(tr, 2, function(x) mean(is.na(x)))
  fzer <- apply(tr, 2, function(x) mean(x == 0, na.rm = TRUE))
  tr[, fna <= thr_na & fzer <= thr_zeros, drop = FALSE]
}
tr_IF <- clean_traits(tr_IF)
tr_SE <- clean_traits(tr_SE)

# T2D status for sample splitting
t2d_IF <- read.delim(file.path(base, "data/traits/Traits_IF.txt"),
                     check.names = FALSE, row.names = 1,
                     na.strings = c("NA", " ", "", "n/a", "N/A"))$T2D
names(t2d_IF) <- rownames(read.delim(file.path(base, "data/traits/Traits_IF.txt"),
                                     check.names = FALSE, row.names = 1))

t2d_SE <- read.delim(file.path(base, "data/traits/Traits_SERUM.txt"),
                     check.names = FALSE, row.names = 1,
                     na.strings = c("NA", " ", "", "n/a", "N/A"))$T2D
names(t2d_SE) <- rownames(read.delim(file.path(base, "data/traits/Traits_SERUM.txt"),
                                     check.names = FALSE, row.names = 1))

IF_cases    <- names(t2d_IF)[!is.na(t2d_IF) & t2d_IF == 1]
IF_controls <- names(t2d_IF)[!is.na(t2d_IF) & t2d_IF == 0]
SE_cases    <- names(t2d_SE)[!is.na(t2d_SE) & t2d_SE == 1]
SE_controls <- names(t2d_SE)[!is.na(t2d_SE) & t2d_SE == 0]

# ==========================================
# CORRELATION HELPER
# ==========================================
compute_cor <- function(proteins, mat, traits, samples) {
  prots_avail <- proteins[proteins %in% rownames(mat)]
  if (length(prots_avail) == 0) return(NULL)

  samp <- intersect(samples, intersect(colnames(mat), rownames(traits)))
  if (length(samp) < 5) return(NULL)

  data   <- t(mat[prots_avail, samp, drop = FALSE])   # samples x proteins
  tr_sub <- traits[samp, , drop = FALSE]

  rho   <- cor(data, tr_sub, use = "pairwise.complete.obs", method = "pearson")
  n_mat <- crossprod(!is.na(data), !is.na(tr_sub))
  t_stat <- rho * sqrt(n_mat - 2) / sqrt(1 - rho^2)
  pval   <- 2 * pt(-abs(t_stat), df = n_mat - 2)
  dimnames(pval) <- dimnames(rho)

  # mask non-significant
  rho_masked <- rho
  rho_masked[pval > pval_thr] <- NA
  list(rho = rho, pval = pval, rho_masked = rho_masked)
}

# ==========================================
# COLOUR PALETTE (Blue-White-Red)
# ==========================================
pal_colors <- colorRampPalette(
  c(rev(brewer.pal(9, "Blues")), "white", brewer.pal(9, "Reds"))
)(200)

# ==========================================
# BUILD LONG DATA FRAME FOR GGPLOT
# ==========================================
# For each group (T2D / control): compute 6 correlation blocks
# Block = protein_group x measured_in_fluid

make_long <- function(group_label, samples_IF, samples_SE) {

  blocks <- list(
    list(prots = if_unique,    mat = mat_IF, tr = tr_IF, samp = samples_IF,
         prot_group = "IF-unique",    fluid_label = "Measured in IF"),
    list(prots = if_unique,    mat = mat_SE, tr = tr_SE, samp = samples_SE,
         prot_group = "IF-unique",    fluid_label = "Measured in SERUM"),
    list(prots = shared,       mat = mat_IF, tr = tr_IF, samp = samples_IF,
         prot_group = "Shared",        fluid_label = "Measured in IF"),
    list(prots = shared,       mat = mat_SE, tr = tr_SE, samp = samples_SE,
         prot_group = "Shared",        fluid_label = "Measured in SERUM"),
    list(prots = serum_unique, mat = mat_SE, tr = tr_SE, samp = samples_SE,
         prot_group = "SERUM-unique",  fluid_label = "Measured in SERUM"),
    list(prots = serum_unique, mat = mat_IF, tr = tr_IF, samp = samples_IF,
         prot_group = "SERUM-unique",  fluid_label = "Measured in IF")
  )

  all_rows <- lapply(blocks, function(b) {
    res <- compute_cor(b$prots, b$mat, b$tr, b$samp)
    if (is.null(res)) return(NULL)
    rho <- res$rho_masked
    df <- data.frame(
      protein    = rep(rownames(rho), times = ncol(rho)),
      trait      = rep(colnames(rho), each  = nrow(rho)),
      r          = as.vector(rho),
      prot_group = b$prot_group,
      fluid      = b$fluid_label,
      group      = group_label,
      stringsAsFactors = FALSE
    )
    df
  })

  do.call(rbind, Filter(Negate(is.null), all_rows))
}

df_T2D  <- make_long("T2D",     IF_cases,    SE_cases)
df_ctrl <- make_long("Control", IF_controls, SE_controls)
df_all  <- rbind(df_T2D, df_ctrl)

# ==========================================
# FIND TRAITS WITH >= 1 SIG CORRELATION
# (across all blocks and both groups)
# ==========================================
sig_traits <- df_all %>%
  filter(!is.na(r)) %>%
  pull(trait) %>%
  unique() %>%
  sort()

df_all <- df_all[df_all$trait %in% sig_traits, ]

# ==========================================
# FACTOR ORDERING
# ==========================================
# Protein order: IF-unique (sorted), then shared, then SERUM-unique
prot_order <- c(
  sort(if_unique[if_unique %in% df_all$protein]),
  sort(shared[shared      %in% df_all$protein]),
  sort(serum_unique[serum_unique %in% df_all$protein])
)
df_all$protein    <- factor(df_all$protein, levels = rev(prot_order))
df_all$trait      <- factor(df_all$trait,   levels = sig_traits)
df_all$fluid      <- factor(df_all$fluid,
                            levels = c("Measured in IF", "Measured in SERUM"))
df_all$prot_group <- factor(df_all$prot_group,
                            levels = c("IF-unique", "Shared", "SERUM-unique"))
df_all$group      <- factor(df_all$group, levels = c("T2D", "Control"))

# ==========================================
# PLOT FUNCTION
# ==========================================
make_heatmap <- function(df_sub, title) {

  # Separator lines between protein groups (at y-breaks)
  group_sizes <- df_sub %>%
    select(protein, prot_group) %>%
    distinct() %>%
    count(prot_group)

  ggplot(df_sub, aes(x = trait, y = protein, fill = r)) +
    geom_tile(color = "grey88", linewidth = 0.25) +
    facet_grid(prot_group ~ fluid, scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colors   = pal_colors,
      limits   = c(-1, 1),
      na.value = "grey96",
      name     = "Pearson r",
      guide    = guide_colorbar(barheight = 10, barwidth = 0.9)
    ) +
    labs(
      title    = title,
      subtitle = paste0("Filled = p \u2264 ", pval_thr,
                        " | white/grey = not significant | rows grouped by DEG origin"),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 7.5),
      axis.text.y      = element_text(size = 7.5),
      strip.text.x     = element_text(face = "bold", size = 9),
      strip.text.y     = element_text(face = "bold", size = 8, angle = 0),
      strip.background = element_rect(fill = "grey92"),
      panel.grid       = element_blank(),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      legend.position  = "right"
    )
}

# ==========================================
# SAVE — one figure per group (T2D / control)
# ==========================================
for (grp in c("T2D", "Control")) {
  df_sub <- df_all[df_all$group == grp, ]
  if (nrow(df_sub) == 0) {
    message("No significant correlations for group: ", grp, " — skipped")
    next
  }
  n_prots  <- length(unique(df_sub$protein))
  n_traits <- length(unique(df_sub$trait))
  h <- max(5, n_prots  * 0.30 + 3)
  w <- max(8, n_traits * 0.45 + 5)

  tag <- paste0("fig6_cross_fluid_trait_cor_", grp)
  p   <- make_heatmap(df_sub,
    title = paste0("DEG \u2013 clinical trait correlations: IF vs SERUM (",
                   grp, " group)"))
  ggsave(file.path(dirOut, paste0(tag, ".pdf")),
         p, width = w, height = h)
  message("Saved: ", tag, ".pdf  (",
          n_prots, " proteins x ", n_traits, " traits)")
}

message("\nDone. Output in: ", dirOut)
