rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# PATHS
# ==========================================
base   <- "~/Documents/projects/ISF/ISF_fede"
dirOut <- file.path(base, "code", "supplementary_analysis", "output")

mat_IF <- read.delim(file.path(base, "data/matrix/IF_joint.txt"),
                     check.names = FALSE, row.names = 1)
tr_IF  <- read.delim(file.path(base, "data/traits/Traits_IF.txt"),
                     check.names = FALSE, row.names = 1,
                     na.strings = c("NA", " ", "", "n/a", "N/A"))

# ==========================================
# DRIVER PROTEINS — IF turquoise
# MM > 0.8 + p.GS < 0.05
# ==========================================
gi <- read.table(
  file.path(base, "code/WGCNA/code/project/ISF/dataset/IF/Results/txtFile/geneInfo.txt"),
  sep = "\t", header = TRUE, quote = "", check.names = FALSE, row.names = 1)
rownames(gi) <- gsub('"', '', rownames(gi))

turq <- gi[gi$moduleColor == "turquoise", ]
drivers <- rownames(turq)[turq$MM.turquoise > 0.8 & turq$p.GS.T2D < 0.05]
cat("Driver proteins:", paste(drivers, collapse = ", "), "\n")

# ==========================================
# LDL-RELATED TRAITS
# ==========================================
ldl_traits <- grep("LDL|ApoB|Chol|TG|chol|aggreg|OX",
                   colnames(tr_IF), value = TRUE, ignore.case = TRUE)
ldl_traits <- ldl_traits[ldl_traits != "T2D"]
cat("LDL traits:", paste(ldl_traits, collapse = ", "), "\n\n")

# ==========================================
# T2D STATUS
# ==========================================
t2d_vec <- tr_IF$T2D
names(t2d_vec) <- rownames(tr_IF)

# ==========================================
# FIND SIGNIFICANT CORRELATIONS IN T2D
# for each driver × LDL trait
# ==========================================
t2d_samp  <- names(t2d_vec)[!is.na(t2d_vec) & t2d_vec == 1]
ctrl_samp <- names(t2d_vec)[!is.na(t2d_vec) & t2d_vec == 0]

sig_pairs <- list()

for (prot in drivers) {
  if (!prot %in% rownames(mat_IF)) next
  for (tr in ldl_traits) {
    if (!tr %in% colnames(tr_IF)) next

    # T2D correlation
    samp <- intersect(t2d_samp, intersect(colnames(mat_IF), rownames(tr_IF)))
    x <- as.numeric(mat_IF[prot, samp])
    y <- as.numeric(tr_IF[samp, tr])
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 5) next

    ct <- cor.test(x[ok], y[ok], method = "pearson")
    if (ct$p.value <= 0.05) {
      sig_pairs[[length(sig_pairs) + 1]] <- list(
        protein = prot, trait = tr,
        r_T2D = ct$estimate, p_T2D = ct$p.value
      )
    }
  }
}

cat("Significant pairs in T2D:", length(sig_pairs), "\n")
if (length(sig_pairs) == 0) {
  message("No significant pairs found — exiting")
  quit(save = "no")
}

# ==========================================
# BUILD LONG DATA FRAME FOR PLOTTING
# ==========================================
all_samp <- intersect(colnames(mat_IF), rownames(tr_IF))

plot_df <- do.call(rbind, lapply(sig_pairs, function(p) {
  x <- as.numeric(mat_IF[p$protein, all_samp])
  y <- as.numeric(tr_IF[all_samp, p$trait])
  g <- ifelse(t2d_vec[all_samp] == 1, "T2D", "Control")
  data.frame(
    protein = p$protein,
    trait   = p$trait,
    expr    = x,
    trait_val = y,
    group   = g,
    r_T2D   = round(p$r_T2D, 2),
    p_T2D   = signif(p$p_T2D, 2),
    stringsAsFactors = FALSE
  )
}))

plot_df <- plot_df[!is.na(plot_df$expr) & !is.na(plot_df$trait_val), ]
plot_df$group <- factor(plot_df$group, levels = c("Control", "T2D"))

# ==========================================
# CORRELATION LABEL PER PANEL
# ==========================================
# Compute r and p for both groups per pair
cor_labels <- plot_df %>%
  group_by(protein, trait, group) %>%
  summarise(
    r = ifelse(n() >= 5,
               round(cor(expr, trait_val, use = "complete.obs"), 2),
               NA_real_),
    p = ifelse(n() >= 5,
               signif(cor.test(expr, trait_val)$p.value, 2),
               NA_real_),
    x_pos = quantile(expr, 0.05, na.rm = TRUE),
    y_pos = max(trait_val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(!is.na(r), paste0("r=", r, ", p=", p), ""),
    sig   = ifelse(!is.na(p) & p <= 0.05, "bold", "plain")
  )

# ==========================================
# COLOURS
# ==========================================
group_cols <- c("Control" = "#1E88E5", "T2D" = "#E53935")

# ==========================================
# PLOT
# ==========================================
make_scatter <- function(prot) {
  df_p  <- plot_df[plot_df$protein == prot, ]
  lbl_p <- cor_labels[cor_labels$protein == prot, ]
  traits_p <- unique(df_p$trait)

  plots <- lapply(traits_p, function(tr) {
    df_t  <- df_p[df_p$trait == tr, ]
    lbl_t <- lbl_p[lbl_p$trait == tr, ]

    # y-label positions: stagger T2D and Control
    y_max <- max(df_t$trait_val, na.rm = TRUE)
    y_rng <- diff(range(df_t$trait_val, na.rm = TRUE))
    lbl_t$y_pos <- c(y_max - y_rng * 0.05,
                     y_max - y_rng * 0.15)[seq_len(nrow(lbl_t))]

    ggplot(df_t, aes(x = expr, y = trait_val, color = group)) +
      geom_point(alpha = 0.75, size = 2) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
      geom_text(data = lbl_t,
                aes(x = x_pos, y = y_pos,
                    label = label, fontface = sig, color = group),
                size = 3, hjust = 0, show.legend = FALSE) +
      scale_color_manual(values = group_cols, name = NULL) +
      labs(
        title = paste0(prot, " \u2014 ", tr),
        x     = paste0(prot, " expression (NPX)"),
        y     = tr
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title      = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
  })
  plots
}

# ==========================================
# SAVE — one PDF per driver protein
# ==========================================
for (prot in drivers) {
  plots_p <- make_scatter(prot)
  if (length(plots_p) == 0) next

  ncols <- min(3, length(plots_p))
  nrows <- ceiling(length(plots_p) / ncols)

  p_combined <- wrap_plots(plots_p, ncol = ncols) +
    plot_annotation(
      title    = paste0("Driver protein: ", prot, " (IF turquoise module)"),
      subtitle = "Scatter: expression vs LDL-related traits | regression lines per group | bold r = p \u2264 0.05",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
    )

  h   <- max(4, nrows * 3.5)
  w   <- ncols * 4
  tag <- paste0("fig7_scatter_", prot)
  ggsave(file.path(dirOut, paste0(tag, ".pdf")),
         p_combined, width = w, height = h)
  message("Saved: ", tag, ".pdf  (", length(plots_p), " panels)")
}

message("\nDone. Output in: ", dirOut)
