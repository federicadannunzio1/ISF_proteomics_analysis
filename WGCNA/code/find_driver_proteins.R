rm(list = ls())

# ==========================================
# PARAMETERS — change before running
# ==========================================
dataset      <- "IF"          # "IF" or "SERUM"
module_color <- "turquoise"   # "turquoise", "brown", etc.
mm_threshold <- 0.8           # 0.8 for IF, 0.7 for SERUM
gs_threshold <- 0.05          # p.GS.T2D significance threshold

# ==========================================
# PATHS
# ==========================================
base     <- "~/Documents/projects/ISF/ISF_fede"
res_path <- file.path(base, "code/WGCNA/code/project/ISF/dataset",
                      dataset, "Results/txtFile/geneInfo.txt")
out_dir  <- file.path(base, "code/WGCNA/code/project/ISF/dataset",
                      dataset, "Results/txtFile")

# ==========================================
# LOAD AND FILTER
# ==========================================
gi     <- read.delim(res_path, check.names = FALSE, row.names = 1)
rownames(gi) <- gsub('"', '', rownames(gi))

module <- gi[gi$moduleColor == module_color, ]
cat(sprintf("Module %s (%s): %d proteins total\n",
            module_color, dataset, nrow(module)))

mm_col <- paste0("MM.", module_color)
if (!mm_col %in% colnames(module)) {
  stop("Column '", mm_col, "' not found in geneInfo. Check module_color.")
}

drivers <- module[
  !is.na(module[[mm_col]])     & module[[mm_col]]     >= mm_threshold &
  !is.na(module$p.GS.T2D)     & module$p.GS.T2D      <= gs_threshold,
]

cat(sprintf("Driver proteins (MM.%s >= %.2f & p.GS.T2D <= %.2f): %d\n",
            module_color, mm_threshold, gs_threshold, nrow(drivers)))

if (nrow(drivers) == 0) {
  message("No driver proteins found with current thresholds.")
} else {
  print(drivers[, c(mm_col, "GS.T2D", "p.GS.T2D")])
}

# ==========================================
# SAVE
# ==========================================
out_file <- file.path(out_dir,
  sprintf("driver_proteins_%s_%s_MM%.2f_GS%.2f.txt",
          dataset, module_color, mm_threshold, gs_threshold))

write.table(drivers, out_file,
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

message("\nSaved: ", out_file,
        "  (", nrow(drivers), " proteins)")
