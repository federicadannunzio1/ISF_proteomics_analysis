rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggvenn)

# ==========================================
# PATHS
# ==========================================
path   <- "~/Documents/projects/ISF/ISF_fede"
dirRes <- file.path(path, "code", "DEGs", "res")
dirOut <- file.path(path, "code", "DEGs", "common")
if (!dir.exists(dirOut)) dir.create(dirOut)

fdr_thr <- 0.05  # soglia DEG

# ==========================================
# Carica risultati DEG e filtra FDR <= 0.05
# ==========================================
deg_if <- read.table(file.path(dirRes, "IF", "DEG.txt"),
                     header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
deg_serum <- read.table(file.path(dirRes, "SERUM", "DEG.txt"),
                        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

sig_if    <- deg_if[deg_if$FDR <= fdr_thr, ]
sig_serum <- deg_serum[deg_serum$FDR <= fdr_thr, ]

message("DEG IF (FDR <= ", fdr_thr, "): ", nrow(sig_if))
message("DEG SERUM (FDR <= ", fdr_thr, "): ", nrow(sig_serum))

# ==========================================
# Proteine comuni
# ==========================================
common_proteins <- intersect(sig_if$protein, sig_serum$protein)
message("Comuni: ", length(common_proteins))

# ==========================================
# Venn diagram
# ==========================================
venn_list <- list(
  IF    = sig_if$protein,
  SERUM = sig_serum$protein
)

pdf(file.path(dirOut, "venn_DEG_IF_vs_SERUM.pdf"), width = 6, height = 5)
print(ggvenn(venn_list,
             fill_color     = c("steelblue", "orange"),
             stroke_color   = c("steelblue", "orange"),
             set_name_color = c("steelblue", "orange"),
             show_percentage = FALSE,
             text_size = 5))
dev.off()
message("Venn salvato.")

# ==========================================
# Tabella riassuntiva
# ==========================================
only_if    <- setdiff(sig_if$protein,    sig_serum$protein)
only_serum <- setdiff(sig_serum$protein, sig_if$protein)

max_len <- max(length(common_proteins), length(only_if), length(only_serum))
pad <- function(v, n) c(v, rep("", n - length(v)))

summary_table <- data.frame(
  common_IF_and_SERUM = pad(common_proteins, max_len),
  only_IF             = pad(only_if,         max_len),
  only_SERUM          = pad(only_serum,      max_len)
)
write.table(summary_table,
            file      = file.path(dirOut, "summary_DEG_IF_vs_SERUM.txt"),
            sep       = "\t", row.names = FALSE, quote = FALSE)

# ==========================================
# Salva lista proteine comuni con dettagli
# ==========================================
common_detail <- merge(
  sig_if[,    c("protein", "FDR", "log2FC", "direction")],
  sig_serum[, c("protein", "FDR", "log2FC", "direction")],
  by     = "protein",
  suffixes = c("_IF", "_SERUM")
)
common_detail <- common_detail[common_detail$protein %in% common_proteins, ]
common_detail <- common_detail[order(common_detail$FDR_IF), ]

write.table(common_detail,
            file      = file.path(dirOut, "common_DEGs.txt"),
            sep       = "\t", row.names = FALSE, quote = FALSE)

message("Tabella comuni salvata: ", nrow(common_detail), " proteine")
message("Output in: ", dirOut)
