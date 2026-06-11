rm(list = ls())
options(stringsAsFactors = FALSE)

library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# PATHS
# ==========================================
path     <- "~/Documents/projects/ISF/ISF_fede"
enr_base <- file.path(path, "code", "WGCNA", "code", "project", "ISF", "dataset")
dirOut   <- file.path(path, "code", "modules_analysis")

# ==========================================
# PARAMETERS
# ==========================================
padj_thr <- 0.05   # significance threshold
top_n    <- 10     # top N terms per module to include

# WGCNA module colors
module_colors <- c(
  turquoise = "turquoise3",
  blue      = "royalblue",
  brown     = "saddlebrown",
  yellow    = "gold2",
  green     = "forestgreen",
  grey      = "grey50"
)

modules_IF    <- c("turquoise", "blue", "brown", "green", "yellow")
modules_SERUM <- c("turquoise", "brown", "yellow")

databases <- c(
  "GO_Biological_Process_2025",
  "GO_Molecular_Function_2025",
  "KEGG_2026",
  "DisGeNET"
)

# ==========================================
# Load enrichment results for one module
# ==========================================
load_enr <- function(fluid, module, db) {
  fname <- paste0(module, "_", db, "_enrichment.xlsx")
  fpath <- file.path(enr_base, fluid, "Results", "functional_enrichment", fname)
  if (!file.exists(fpath)) return(NULL)
  df <- suppressMessages(read_excel(fpath))
  df <- df[!is.na(df$Adjusted.P.value) & df$Adjusted.P.value <= padj_thr, ]
  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$Adjusted.P.value), ]
  df <- head(df, top_n)
  data.frame(
    Term   = df$Term,
    padj   = df$Adjusted.P.value,
    score  = df$Combined.Score,
    fluid  = fluid,
    module = module,
    stringsAsFactors = FALSE
  )
}

# ==========================================
# Main loop: one plot per database
# ==========================================
for (db in databases) {

  message("Processing: ", db)

  chunks <- list()
  for (mod in modules_IF)    chunks[[length(chunks)+1]] <- load_enr("IF",    mod, db)
  for (mod in modules_SERUM) chunks[[length(chunks)+1]] <- load_enr("SERUM", mod, db)

  df_all <- do.call(rbind, Filter(Negate(is.null), chunks))

  if (is.null(df_all) || nrow(df_all) == 0) {
    message("  No significant terms for ", db, " — skipped")
    next
  }

  # Truncate long term names
  df_all$Term_short <- ifelse(
    nchar(df_all$Term) > 55,
    paste0(substr(df_all$Term, 1, 52), "..."),
    df_all$Term
  )

  # ---- Shared Y axis ----
  # Order terms: those enriched in BOTH fluids first, then by min padj
  term_fluids <- df_all %>%
    group_by(Term_short) %>%
    summarise(
      n_fluids = n_distinct(fluid),
      min_padj = min(padj),
      .groups  = "drop"
    ) %>%
    arrange(desc(n_fluids), min_padj)

  shared_terms  <- term_fluids$Term_short[term_fluids$n_fluids == 2]
  unique_terms  <- term_fluids$Term_short[term_fluids$n_fluids == 1]
  term_order    <- c(shared_terms, unique_terms)

  df_all$Term_short <- factor(df_all$Term_short, levels = rev(term_order))

  present_colors <- module_colors[names(module_colors) %in% df_all$module]

  # ---- Horizontal barplot with shared Y axis ----
  # Facet by fluid, free_x (each panel its own x range),
  # fixed Y so both panels show the same terms — empty gap where absent
  p <- ggplot(df_all, aes(x = -log10(padj), y = Term_short, fill = module)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = -log10(padj_thr), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    scale_fill_manual(values = present_colors, name = "WGCNA module") +
    facet_wrap(~ fluid, ncol = 2, scales = "free_x") +
    labs(
      title    = paste0("Functional enrichment comparison: IF vs SERUM — ", db),
      subtitle = paste0("Top ", top_n, " terms per module | adj. p <= ", padj_thr,
                        " | Terms enriched in both fluids shown at top"),
      x = "-log10(adj. p-value)", y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text         = element_text(face = "bold", size = 13),
      axis.text.y        = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      legend.position    = "right",
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 9)
    )

  n_terms <- length(levels(df_all$Term_short))
  h <- max(5, n_terms * 0.25)

  ggsave(file.path(dirOut, paste0("enrichment_comparison_", db, ".pdf")),
         plot = p, width = 12, height = h)
  ggsave(file.path(dirOut, paste0("enrichment_comparison_", db, ".png")),
         plot = p, width = 12, height = h, dpi = 300)

  n_shared <- length(shared_terms)
  message("  Saved: enrichment_comparison_", db, ".pdf  (",
          n_terms, " terms, ", n_shared, " shared between IF and SERUM)")
}

message("\nDone. Output in: ", dirOut)
