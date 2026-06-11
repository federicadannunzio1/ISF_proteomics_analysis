rm(list = ls())
options(stringsAsFactors = FALSE)

library(WGCNA)
library(ggplot2)
library(ggrepel)

enableWGCNAThreads()

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
base_path  <- "~/Documents/projects/ISF/ISF_fede"
wgcna_base <- file.path(base_path, "code", "WGCNA", "code", "project", "ISF", "dataset")

path_net_IF   <- file.path(wgcna_base, "IF",    "Results", "RData", "networkConstruction.RData")
path_net_SE   <- file.path(wgcna_base, "SERUM", "Results", "RData", "networkConstruction.RData")
path_data_IF  <- file.path(wgcna_base, "IF",    "Results", "RData", "dataInput.RData")
path_data_SE  <- file.path(wgcna_base, "SERUM", "Results", "RData", "dataInput.RData")

out_dir <- file.path(base_path, "code", "WGCNA", "res", "module_preservation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────
n_perms    <- 200   # number of permutations (200 standard; use 50 for quick test)
rand_seed  <- 42
z_thr_low  <- 2     # Z < 2: not preserved
z_thr_high <- 10    # Z > 10: highly preserved

# ─────────────────────────────────────────────────────────────────────────────
# Load WGCNA results
# ─────────────────────────────────────────────────────────────────────────────
# IF
load(path_data_IF)
datExpr_IF    <- datExpr
rm(datExpr, datTraits)

load(path_net_IF)
moduleColors_IF <- moduleColors
rm(geneTree, MEs, moduleColors, moduleLabels)

# SERUM
load(path_data_SE)
datExpr_SE    <- datExpr
rm(datExpr, datTraits)

load(path_net_SE)
moduleColors_SE <- moduleColors
rm(geneTree, MEs, moduleColors, moduleLabels)

message("IF:    ", nrow(datExpr_IF), " samples x ", ncol(datExpr_IF), " proteins")
message("SERUM: ", nrow(datExpr_SE), " samples x ", ncol(datExpr_SE), " proteins")

# Proteins in common between the two panels
common_prots <- intersect(colnames(datExpr_IF), colnames(datExpr_SE))
message("Proteins in common between IF and SERUM panels: ", length(common_prots))
# Note: modulePreservation() uses only common proteins automatically

# ─────────────────────────────────────────────────────────────────────────────
# Run modulePreservation: IF as reference, SERUM as test
# ─────────────────────────────────────────────────────────────────────────────
# multiData: list of expression matrices (samples x proteins), one per dataset
# multiColor: list of module color vectors, matching columns of multiData
# The function tests whether the co-expression structure of each IF module
# is reproduced in the SERUM dataset using permutation-based Z-scores.

message("\n--- Module preservation: IF (reference) -> SERUM (test) ---")

multiExpr <- list(
  IF    = list(data = datExpr_IF),
  SERUM = list(data = datExpr_SE)
)
multiColor <- list(
  IF    = moduleColors_IF,
  SERUM = moduleColors_SE
)

set.seed(rand_seed)
mp_IF_ref <- modulePreservation(
  multiData          = multiExpr,
  multiColor         = multiColor,
  dataIsExpr         = TRUE,
  referenceNetworks  = 1,          # IF is the reference
  testNetworks       = 2,          # SERUM is the test
  nPermutations      = n_perms,
  randomSeed         = rand_seed,
  quickCor           = 0,
  verbose            = 3
)

# ─────────────────────────────────────────────────────────────────────────────
# Run modulePreservation: SERUM as reference, IF as test
# ─────────────────────────────────────────────────────────────────────────────
message("\n--- Module preservation: SERUM (reference) -> IF (test) ---")

set.seed(rand_seed)
mp_SE_ref <- modulePreservation(
  multiData          = multiExpr,
  multiColor         = multiColor,
  dataIsExpr         = TRUE,
  referenceNetworks  = 2,          # SERUM is the reference
  testNetworks       = 1,          # IF is the test
  nPermutations      = n_perms,
  randomSeed         = rand_seed,
  quickCor           = 0,
  verbose            = 3
)

# Save RData for downstream use
save(mp_IF_ref, mp_SE_ref, file = file.path(out_dir, "modulePreservation.RData"))
message("RData saved.")

# ─────────────────────────────────────────────────────────────────────────────
# Extract preservation statistics
# ─────────────────────────────────────────────────────────────────────────────
extract_preservation <- function(mp, ref_name, test_name) {
  # Zsummary: overall preservation Z-score (summary of density + connectivity)
  # medianRank: rank-based preservation (lower = better preserved)
  ref_idx  <- paste0("ref.", ref_name)
  test_idx <- paste0("inColumnsAlsoPresentIn.", test_name)

  Z          <- mp$preservation$Z[[ref_idx]][[test_idx]]
  median_rk  <- mp$preservation$observed[[ref_idx]][[test_idx]]

  modules    <- rownames(Z)
  mod_size   <- mp$preservation$Z[[ref_idx]][[test_idx]][, "moduleSize"]

  df <- data.frame(
    module      = modules,
    module_size = as.integer(mod_size),
    Zsummary    = Z[, "Zsummary.pres"],
    medianRank  = median_rk[, "medianRank.pres"],
    reference   = ref_name,
    test        = test_name,
    stringsAsFactors = FALSE
  )

  # Exclude gold (all proteins pooled) and grey (unassigned) from interpretation
  df <- df[!df$module %in% c("gold", "grey"), ]
  df[order(df$Zsummary), ]
}

res_IF_ref <- extract_preservation(mp_IF_ref, ref_name = "IF",    test_name = "SERUM")
res_SE_ref <- extract_preservation(mp_SE_ref, ref_name = "SERUM", test_name = "IF")

# Label preservation category
label_pres <- function(z) {
  ifelse(z < z_thr_low, "not preserved",
         ifelse(z < z_thr_high, "moderately preserved", "highly preserved"))
}
res_IF_ref$preservation <- label_pres(res_IF_ref$Zsummary)
res_SE_ref$preservation <- label_pres(res_SE_ref$Zsummary)

# Save tables
write.table(res_IF_ref,
            file = file.path(out_dir, "preservation_IF_ref.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(res_SE_ref,
            file = file.path(out_dir, "preservation_SERUM_ref.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\n=== IF modules preservation in SERUM ===")
print(res_IF_ref[, c("module", "module_size", "Zsummary", "preservation")])

message("\n=== SERUM modules preservation in IF ===")
print(res_SE_ref[, c("module", "module_size", "Zsummary", "preservation")])

# ─────────────────────────────────────────────────────────────────────────────
# Plot: Zsummary vs module size
# Standard WGCNA preservation plot, one panel per direction
# ─────────────────────────────────────────────────────────────────────────────
make_preservation_plot <- function(df, title) {

  # Use actual module color as point color where possible
  # (module names in WGCNA are color names)
  valid_colors <- tryCatch(
    { col2rgb(df$module); df$module },
    error = function(e) rep("grey50", nrow(df))
  )

  ggplot(df, aes(x = module_size, y = Zsummary, label = module)) +
    geom_hline(yintercept = z_thr_low,  linetype = "dashed", color = "red",    linewidth = 0.8) +
    geom_hline(yintercept = z_thr_high, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
    geom_point(aes(color = module), size = 5) +
    scale_color_manual(
      values = setNames(valid_colors, df$module),
      guide  = "none"
    ) +
    geom_text_repel(size = 4, fontface = "bold",
                    box.padding = 0.5, point.padding = 0.3,
                    max.overlaps = Inf) +
    annotate("text", x = Inf, y = z_thr_low  + 0.3, hjust = 1.1,
             label = "Z = 2 (not preserved)", color = "red", size = 3.2) +
    annotate("text", x = Inf, y = z_thr_high + 0.3, hjust = 1.1,
             label = "Z = 10 (highly preserved)", color = "darkgreen", size = 3.2) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    labs(title    = title,
         subtitle = paste0("n_permutations = ", n_perms,
                           " | grey and gold modules excluded"),
         x = "Module size (number of proteins)",
         y = "Zsummary preservation") +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))
}

p_IF_ref <- make_preservation_plot(
  res_IF_ref,
  title = "Module preservation: IF modules in SERUM"
)
p_SE_ref <- make_preservation_plot(
  res_SE_ref,
  title = "Module preservation: SERUM modules in IF"
)

ggsave(file.path(out_dir, "preservation_IF_ref.pdf"),
       plot = p_IF_ref, width = 7, height = 5)
ggsave(file.path(out_dir, "preservation_SERUM_ref.pdf"),
       plot = p_SE_ref, width = 7, height = 5)

# ─────────────────────────────────────────────────────────────────────────────
# Combined plot: both directions side by side
# ─────────────────────────────────────────────────────────────────────────────
library(patchwork)
p_combined <- (p_IF_ref | p_SE_ref) +
  plot_annotation(
    title    = "WGCNA Module Preservation: IF <-> SERUM",
    subtitle = "Each point = one module. Dashed lines: Z=2 (not preserved) and Z=10 (highly preserved).",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

ggsave(file.path(out_dir, "preservation_combined.pdf"),
       plot = p_combined, width = 13, height = 5)
ggsave(file.path(out_dir, "preservation_combined.png"),
       plot = p_combined, width = 13, height = 5, dpi = 300)

message("\nAll output saved in: ", out_dir)
message("Done.")
