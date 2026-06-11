rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# PATHS
# ==========================================
path    <- "~/Documents/projects/ISF/ISF_fede"
dirOut  <- file.path(path, "code", "supplementary_analysis", "output")
wgcna   <- file.path(path, "code", "WGCNA", "code", "project", "ISF", "dataset")
deg_dir <- file.path(path, "code", "DEGs", "res")

module_colors <- c(
  turquoise = "turquoise3", blue = "royalblue", brown = "saddlebrown",
  yellow = "gold2", green = "forestgreen", grey = "grey70"
)

# ==========================================
# LOAD DATA
# ==========================================
load_geneinfo <- function(fluid) {
  gi <- read.table(file.path(wgcna, fluid, "Results", "txtFile", "geneInfo.txt"),
                   sep = "\t", header = TRUE, quote = "", check.names = FALSE,
                   row.names = 1)
  rownames(gi) <- gsub('"', '', rownames(gi))
  gi$protein <- rownames(gi)
  gi$fluid   <- fluid
  gi
}

gi_IF    <- load_geneinfo("IF")
gi_SERUM <- load_geneinfo("SERUM")

deg_IF    <- read.delim(file.path(deg_dir, "IF",    "DEG_res", "DEG.txt"), sep = "\t")
deg_SERUM <- read.delim(file.path(deg_dir, "SERUM", "DEG_res", "DEG.txt"), sep = "\t")

if_unique    <- deg_IF$protein[!(deg_IF$protein    %in% deg_SERUM$protein)]
serum_unique <- deg_SERUM$protein[!(deg_SERUM$protein %in% deg_IF$protein)]

# ==========================================
# FIGURE A — kME boxplot per module, IF vs SERUM
# (own-module MM only, grey excluded)
# ==========================================
extract_own_kme <- function(gi) {
  mods <- setdiff(unique(gi$moduleColor), "grey")
  do.call(rbind, lapply(mods, function(m) {
    members <- gi[gi$moduleColor == m, ]
    mm_col  <- paste0("MM.", m)
    if (!mm_col %in% names(members)) return(NULL)
    data.frame(
      protein = members$protein,
      module  = m,
      fluid   = members$fluid[1],
      kME     = members[[mm_col]],
      stringsAsFactors = FALSE
    )
  }))
}

kme_IF    <- extract_own_kme(gi_IF)
kme_SERUM <- extract_own_kme(gi_SERUM)
kme_all   <- rbind(kme_IF, kme_SERUM)

# Order modules: IF first, then SERUM; within each fluid by median kME desc
mod_order <- kme_all %>%
  group_by(fluid, module) %>%
  summarise(med = median(kME), .groups = "drop") %>%
  arrange(fluid, desc(med)) %>%
  mutate(label = paste0(module, "\n(", fluid, ")"))

kme_all <- kme_all %>%
  left_join(mod_order %>% select(fluid, module, label), by = c("fluid", "module"))

kme_all$label  <- factor(kme_all$label,  levels = mod_order$label)
kme_all$fluid  <- factor(kme_all$fluid,  levels = c("IF", "SERUM"))

fill_vals <- setNames(
  module_colors[mod_order$module],
  mod_order$label
)

p_kme <- ggplot(kme_all, aes(x = label, y = kME, fill = label)) +
  geom_boxplot(outlier.size = 1.5, width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40",
             linewidth = 0.5) +
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "grey60",
             linewidth = 0.5) +
  facet_grid(. ~ fluid, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fill_vals, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = Inf, y = 0.81, label = "MM = 0.8",
           hjust = 1.1, size = 3, color = "grey40") +
  labs(
    title    = "Intramodular connectivity (kME) — IF vs SERUM",
    subtitle = "Each box = distribution of own-module MM for all module members | grey excluded",
    x = NULL, y = "Module Membership (kME)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(size = 9),
    strip.text      = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92"),
    panel.grid.major.x = element_blank(),
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

# ==========================================
# FIGURE B — DEG × module overlap: stacked bar
# ==========================================
# IF
if_mod <- data.frame(
  protein = if_unique,
  module  = ifelse(if_unique %in% rownames(gi_IF),
                   gi_IF[if_unique[if_unique %in% rownames(gi_IF)], "moduleColor"],
                   "not detected"),
  fluid   = "IF",
  stringsAsFactors = FALSE
)
# handle proteins not in WGCNA
if_mod$module[!(if_unique %in% rownames(gi_IF))] <- "not detected"

# SERUM
serum_mod <- data.frame(
  protein = serum_unique,
  module  = ifelse(serum_unique %in% rownames(gi_SERUM),
                   gi_SERUM[serum_unique[serum_unique %in% rownames(gi_SERUM)], "moduleColor"],
                   "not detected"),
  fluid   = "SERUM",
  stringsAsFactors = FALSE
)
serum_mod$module[!(serum_unique %in% rownames(gi_SERUM))] <- "not detected"

deg_mod <- rbind(if_mod, serum_mod)

# Count per fluid × module
deg_counts <- deg_mod %>%
  count(fluid, module) %>%
  group_by(fluid) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# Module order for fill
all_mods <- c("turquoise","blue","brown","yellow","green","grey","not detected")
deg_counts$module <- factor(deg_counts$module, levels = rev(all_mods))

fill_deg <- c(module_colors, "not detected" = "#BDBDBD")
fill_deg <- fill_deg[levels(deg_counts$module)]

p_deg <- ggplot(deg_counts, aes(x = fluid, y = n, fill = module)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = fill_deg, name = "WGCNA module",
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = seq(0, 20, 2)) +
  labs(
    title    = "Fluid-unique DEGs — WGCNA module distribution",
    subtitle = "n = proteins per module | grey = unassigned",
    x = NULL, y = "Number of DEGs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    legend.position    = "right"
  )

# ==========================================
# SAVE
# ==========================================
ggsave(file.path(dirOut, "fig5a_kME_distribution_IF_vs_SERUM.pdf"),
       p_kme, width = 10, height = 5)
message("Saved: fig5a_kME_distribution_IF_vs_SERUM.pdf")

ggsave(file.path(dirOut, "fig5b_DEG_module_overlap.pdf"),
       p_deg, width = 6, height = 5)
message("Saved: fig5b_DEG_module_overlap.pdf")

p_full <- p_kme / p_deg +
  plot_annotation(
    title = "Network architecture: IF vs SERUM",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )
ggsave(file.path(dirOut, "fig5_network_architecture.pdf"),
       p_full, width = 11, height = 10)
message("Saved: fig5_network_architecture.pdf")

# ==========================================
# PRINT SUMMARY
# ==========================================
cat("\n=== DEG-module overlap summary ===\n")
for (fl in c("IF","SERUM")) {
  cat(sprintf("\n%s-unique DEGs (%d total):\n", fl,
              sum(deg_counts$fluid == fl & deg_counts$module != "not detected",
                  na.rm = TRUE) |> {\(x) deg_counts$n[deg_counts$fluid==fl] |> sum()}()))
  sub <- deg_counts[deg_counts$fluid == fl, ]
  sub <- sub[order(-sub$n), ]
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-15s %d proteins (%.0f%%)\n",
                as.character(sub$module[i]), sub$n[i], sub$pct[i]))
  }
}
