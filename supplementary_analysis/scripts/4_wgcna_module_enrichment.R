rm(list = ls())
options(stringsAsFactors = FALSE)

library(readxl)
library(ggplot2)
library(dplyr)

# ==========================================
# PATHS & PARAMETERS
# ==========================================
path     <- "~/Documents/projects/ISF/ISF_fede"
wgcna_base <- file.path(path, "code", "WGCNA", "code", "project", "ISF", "dataset")
dirOut   <- file.path(path, "code", "supplementary_analysis", "output")
dir.create(dirOut, showWarnings = FALSE)

padj_thr       <- 0.05
top_n_per_mod  <- 10   # top terms per module before union
top_n_global   <- 25   # max terms shown in plot A

databases <- c(
  "GO_Biological_Process_2025",
  "GO_Molecular_Function_2025",
  "KEGG_2026",
  "DisGeNET"
)

modules_IF   <- c("turquoise", "blue", "brown", "yellow", "green")
modules_SERUM <- c("turquoise", "brown", "yellow")

# ggplot-friendly colours for module labels / fill
module_colors <- c(
  turquoise = "turquoise3",
  blue      = "royalblue",
  brown     = "saddlebrown",
  yellow    = "gold2",
  green     = "forestgreen"
)

# ==========================================
# CARDIAC EXCLUSION
# ==========================================
cardiac_exclude_kw <- c(
  "cardiac atrium", "cardiac ventricle", "ventricular cardiac",
  "cardiac ventricle morpho", "cardiac ventricle dev",
  "ventricular septum", "atrial", "atrioventricular",
  "cardiac muscle", "myocardium", "ventricular trabecula", "heart trabecula",
  "trabecula",
  "heart morphogenesis", "heart looping", "heart tube", "heart left.right",
  "determination of heart", "embryonic heart", "cardiac cell",
  "aortic valve", "pulmonary valve", "mitral valve", "outflow tract",
  "endocardial cushion", "endocardial", "cardiac epithelial",
  "cardiac epithelial to mesenchymal",
  "pericardium", "cardioblast",
  "myoblast differentiation"
)

is_cardiac <- function(term) {
  t <- tolower(term)
  any(sapply(cardiac_exclude_kw, function(k) grepl(k, t, fixed = TRUE)))
}

# ==========================================
# LOAD HELPER
# ==========================================
load_module_enr <- function(fluid, module, db) {
  fpath <- file.path(wgcna_base, fluid, "Results", "functional_enrichment",
                     paste0(module, "_", db, "_enrichment.xlsx"))
  if (!file.exists(fpath)) return(NULL)
  df <- suppressMessages(read_excel(fpath))
  df <- df[!is.na(df$Adjusted.P.value) & df$Adjusted.P.value <= padj_thr, ]
  if (nrow(df) == 0) return(NULL)
  df <- df[!sapply(df$Term, is_cardiac), ]
  if (nrow(df) == 0) return(NULL)
  df$module <- module
  df$fluid  <- fluid
  df
}

# ==========================================
# BIOLOGICAL CATEGORY (for Plot B colour)
# ==========================================
categorize_term <- function(term) {
  t <- tolower(term)

  cardiac_kw <- c("cardiac atrium", "cardiac ventricle", "cardiac muscle tissue",
                  "cardiac epithelial", "ventricular trabecula", "heart trabecula",
                  "endocardial cushion", "heart looping", "heart tube",
                  "heart morphogenesis", "ventricular cardiac", "myocardium morphogenesis",
                  "atrioventricular", "myoblast differentiation", "heart left.right",
                  "embryonic heart")
  if (any(sapply(cardiac_kw, function(k) grepl(k, t, fixed = TRUE))))
    return("Cardiac Development")

  vasc_kw <- c("vascular", "vessel", "artery", "arteri", "endotheli", "smooth muscle",
               "vasculo", "angiogen", "mesenchyme", "extracellular matrix",
               "bmp signal", "tgf", "fibrosis", "collagen", "fibronectin",
               "epithelial to mesenchymal", "cell adhesion", "ecm")
  if (any(sapply(vasc_kw, function(k) grepl(k, t, fixed = TRUE))))
    return("Vascular/ECM")

  imm_kw <- c("inflamm", "cytokine", "interleukin", "leukocyte", "immune",
              "chemokine", "neuroinflammatory", "glial", "gliogen",
              "interferon", "lymphocyte", "monocyte", "macrophage",
              "neutrophil", "toll-like", "innate", "adaptive", "t cell",
              "b cell", "nk cell", "mast cell", "wounding", "acute")
  if (any(sapply(imm_kw, function(k) grepl(k, t, fixed = TRUE))))
    return("Immune/Inflammatory")

  met_kw <- c("lipid", "fat ", "adipogen", "cholesterol", "insulin",
              "glucose", "glucocorticoid", "corticosteroid", "steroid",
              "metabol", "lipoprotein", "fatty acid", "triglycerid",
              "fat cell", "lipid storage", "lipid localiz")
  if (any(sapply(met_kw, function(k) grepl(k, t, fixed = TRUE))))
    return("Metabolic/Lipid")

  return("Other/Signaling")
}

cat_colors <- c(
  "Cardiac Development" = "#B0BEC5",
  "Vascular/ECM"        = "#E64A19",
  "Immune/Inflammatory" = "#7B1FA2",
  "Metabolic/Lipid"     = "#388E3C",
  "Other/Signaling"     = "#1976D2"
)

# short term label helper
shorten <- function(x, n = 55)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)

# ==========================================
# PLOT A — multi-module dotplot
# rows = terms, columns = modules
# size = -log10(padj), colour = Combined.Score
# ==========================================
make_dotplot <- function(fluid, db) {

  mods <- if (fluid == "IF") modules_IF else modules_SERUM

  # Load and combine
  dfs <- lapply(mods, function(m) load_module_enr(fluid, m, db))
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(NULL)
  df_all <- do.call(rbind, dfs)

  # Select union of top N per module
  top_terms <- df_all %>%
    group_by(module) %>%
    slice_min(Adjusted.P.value, n = top_n_per_mod, with_ties = FALSE) %>%
    ungroup() %>%
    pull(Term) %>%
    unique()

  # From the union, keep top_n_global by best padj across modules
  best_padj <- df_all %>%
    filter(Term %in% top_terms) %>%
    group_by(Term) %>%
    summarise(best_padj = min(Adjusted.P.value), .groups = "drop") %>%
    arrange(best_padj) %>%
    slice_head(n = top_n_global)

  selected_terms <- best_padj$Term

  # Build full grid: all selected_terms × all modules
  grid <- expand.grid(Term = selected_terms, module = mods,
                      stringsAsFactors = FALSE)
  df_plot <- merge(grid, df_all[, c("Term", "module", "Adjusted.P.value", "Combined.Score")],
                   by = c("Term", "module"), all.x = TRUE)

  # Short term labels
  df_plot$Term_short <- shorten(df_plot$Term)

  # Order terms: most significant at top
  term_order <- best_padj %>%
    mutate(Term_short = shorten(Term)) %>%
    pull(Term_short) %>%
    rev()  # rev so most significant ends up at top of y axis

  df_plot$Term_short <- factor(df_plot$Term_short, levels = term_order)
  df_plot$module     <- factor(df_plot$module, levels = mods)

  df_plot$log10p <- -log10(df_plot$Adjusted.P.value)

  # Module label colours for x axis
  x_colors <- module_colors[levels(df_plot$module)]

  p <- ggplot(df_plot[!is.na(df_plot$log10p), ],
              aes(x = module, y = Term_short,
                  size = log10p, colour = Combined.Score)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(
      range  = c(2, 9),
      name   = "-log10(adj. p)",
      breaks = c(2, 5, 10, 20)
    ) +
    scale_colour_gradient(
      low  = "#FFF9C4",
      high = "#E53935",
      name = "Combined\nScore"
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title    = paste0(fluid, " — ", gsub("_", " ", db)),
      subtitle = paste0("Top ", top_n_global, " terms (union of top ", top_n_per_mod,
                        " per module) | adj. p \u2264 ", padj_thr,
                        " | cardiac terms excluded"),
      x = "WGCNA module", y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(colour = x_colors, face = "bold", size = 11),
      axis.text.y      = element_text(size = 8),
      panel.grid.major = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      legend.position  = "right"
    )

  n_terms <- length(selected_terms)
  n_mods  <- length(mods)
  list(plot = p, n_terms = n_terms, n_mods = n_mods)
}

# ==========================================
# PLOT B — two-sided: any module pair across fluids
# fluid1/mod1 shown on left (negative x)
# fluid2/mod2 shown on right (positive x)
# ==========================================
make_twosided_modules <- function(db,
                                  fluid1 = "IF",   mod1 = "turquoise",
                                  fluid2 = "SERUM", mod2 = "brown") {

  df1 <- load_module_enr(fluid1, mod1, db)
  df2 <- load_module_enr(fluid2, mod2, db)

  # Use a single "side" label to distinguish the two groups
  if (!is.null(df1)) df1$side <- paste0(fluid1, " ", mod1)
  if (!is.null(df2)) df2$side <- paste0(fluid2, " ", mod2)

  df_all <- do.call(rbind, Filter(Negate(is.null), list(df1, df2)))
  if (is.null(df_all) || nrow(df_all) == 0) return(NULL)

  label_left  <- paste0(fluid1, " ", mod1)
  label_right <- paste0(fluid2, " ", mod2)
  color_left  <- module_colors[mod1]
  color_right <- module_colors[mod2]

  # Short labels + deduplicate
  df_all$Term_short <- shorten(df_all$Term, 60)
  df_all <- df_all %>%
    group_by(Term_short, side) %>%
    slice_min(Adjusted.P.value, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Top 20 per side
  top_terms <- df_all %>%
    group_by(side) %>%
    slice_min(Adjusted.P.value, n = 20, with_ties = FALSE) %>%
    ungroup() %>%
    pull(Term_short) %>%
    unique()

  df_all <- df_all[df_all$Term_short %in% top_terms, ]

  # Signed x: left side = negative, right side = positive
  df_all$x_val <- ifelse(
    df_all$side == label_left,
    -(-log10(df_all$Adjusted.P.value)),
     -log10(df_all$Adjusted.P.value)
  )

  # Assign categories
  df_all$category <- sapply(df_all$Term_short, categorize_term)

  # Term order: shared first, then left-only, then right-only; within each by |x| desc
  term_info <- df_all %>%
    group_by(Term_short) %>%
    summarise(
      n_sides   = n_distinct(side),
      max_score = max(abs(x_val)),
      .groups   = "drop"
    ) %>%
    arrange(desc(n_sides), desc(max_score))

  df_all$Term_short <- factor(df_all$Term_short,
                              levels = rev(term_info$Term_short))

  x_lim <- max(abs(df_all$x_val)) * 1.1
  n_terms <- nlevels(df_all$Term_short)

  p <- ggplot(df_all, aes(x = x_val, y = Term_short, fill = category)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "grey20", linewidth = 0.6) +
    scale_fill_manual(values = cat_colors, name = "Category", drop = TRUE) +
    scale_x_continuous(
      limits = c(-x_lim, x_lim),
      labels = function(x) round(abs(x), 1),
      name   = "-log10(adj. p-value)"
    ) +
    annotate("text", x = -x_lim * 0.95, y = Inf,
             label = label_left, fontface = "bold", color = color_left,
             hjust = 0, vjust = 1.5, size = 4.5) +
    annotate("text", x =  x_lim * 0.95, y = Inf,
             label = label_right, fontface = "bold", color = color_right,
             hjust = 1, vjust = 1.5, size = 4.5) +
    labs(
      title    = paste0(label_left, " vs ", label_right, " — ", gsub("_", " ", db)),
      subtitle = paste0("Top 20 per module | adj. p \u2264 ", padj_thr,
                        " | cardiac terms excluded"),
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y        = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 9, color = "grey40"),
      legend.position    = "right"
    )

  list(plot = p, n_terms = n_terms)
}

# ==========================================
# MAIN LOOP
# ==========================================
for (db in databases) {
  message("\n========== ", db, " ==========")

  # --- Plot A: IF ---
  res <- make_dotplot("IF", db)
  if (!is.null(res)) {
    h <- max(5, res$n_terms * 0.35)
    w <- max(6, res$n_mods  * 1.2 + 4)
    tag <- paste0("fig4a_dotplot_IF_", db)
    ggsave(file.path(dirOut, paste0(tag, ".pdf")),
           res$plot, width = w, height = h)
    message("  [IF dotplot]    Saved: ", tag,
            "  (", res$n_terms, " terms, ", res$n_mods, " modules)")
  } else {
    message("  [IF dotplot]    No data — skipped")
  }

  # --- Plot A: SERUM ---
  res <- make_dotplot("SERUM", db)
  if (!is.null(res)) {
    h <- max(5, res$n_terms * 0.35)
    w <- max(6, res$n_mods  * 1.2 + 4)
    tag <- paste0("fig4a_dotplot_SERUM_", db)
    ggsave(file.path(dirOut, paste0(tag, ".pdf")),
           res$plot, width = w, height = h)
    message("  [SERUM dotplot] Saved: ", tag,
            "  (", res$n_terms, " terms, ", res$n_mods, " modules)")
  } else {
    message("  [SERUM dotplot] No data — skipped")
  }

  # --- Plot B: IF turquoise vs SERUM brown ---
  res <- make_twosided_modules(db,
                               fluid1 = "IF",    mod1 = "turquoise",
                               fluid2 = "SERUM", mod2 = "brown")
  if (!is.null(res)) {
    h <- max(4, res$n_terms * 0.32)
    tag <- paste0("fig4b_IF_turquoise_vs_SERUM_brown_", db)
    ggsave(file.path(dirOut, paste0(tag, ".pdf")),
           res$plot, width = 12, height = h)
    message("  [IF turq vs SERUM brown] Saved: ", tag,
            "  (", res$n_terms, " terms)")
  } else {
    message("  [IF turq vs SERUM brown] No data — skipped")
  }
}

message("\nDone. Output in: ", dirOut)
