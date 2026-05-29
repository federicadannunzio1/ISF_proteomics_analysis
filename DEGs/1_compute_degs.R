rm(list = ls())
options(stringsAsFactors = FALSE)

library(pheatmap)
library(ggplot2)
library(ggrepel)

# ==========================================
# PATHS
# ==========================================
path        <- "~/Documents/projects/ISF/ISF_fede"
path_matrix <- file.path(path, "data", "matrix")
path_list   <- file.path(path, "data", "list")
dirRes      <- file.path(path, "code", "DEGs", "res")

# ==========================================
# PARAMETERS
# ==========================================
prc              <- 0.1   # remove bottom prc fraction of proteins by IQR (low-variation filter)
soglia_log2fc    <- 0     # log2FC filter for analysis: 0 = no filter, all proteins kept
soglia_pval      <- 0.05  # FDR-adjusted p-value threshold
soglia_log2fc_plot <- 0.2 # log2FC threshold shown as blue lines in volcano (only visual)

# ==========================================
# ANALYSIS LOOP: IF vs Control, SERUM vs Control
# ==========================================

for (fluid in c("IF", "SERUM")) {

  message("\n========== DEG analysis: T2D vs Control in ", fluid, " ==========")

  # ------------------------------------------------------------------
  # 1. Build joint sample lists: union of CMET and INF lists per fluid
  # ------------------------------------------------------------------
  # File naming: CMET_IF_list_case.txt / INF_IF_list_case.txt (for IF)
  #              CMET_serum_list_case.txt / INF_serum_list_case.txt (for SERUM)
  fluid_tag <- ifelse(fluid == "SERUM", "serum", "IF")

  case_files <- list.files(path_list,
                           pattern    = paste0(".*_", fluid_tag, "_list_case\\.txt$"),
                           full.names = TRUE)
  ctrl_files <- list.files(path_list,
                           pattern    = paste0(".*_", fluid_tag, "_list_ctrl\\.txt$"),
                           full.names = TRUE)

  stopifnot("No case list files found" = length(case_files) >= 1,
            "No ctrl list files found" = length(ctrl_files) >= 1)

  case_ids_joint <- unique(trimws(unlist(lapply(case_files, readLines))))
  ctrl_ids_joint <- unique(trimws(unlist(lapply(ctrl_files, readLines))))

  # Save joint lists to disk (union of CMET + INF panels)
  write.table(case_ids_joint,
              file      = file.path(path_list, paste0(fluid, "_joint_list_case.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(ctrl_ids_joint,
              file      = file.path(path_list, paste0(fluid, "_joint_list_ctrl.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  message("Joint ", fluid, " case list saved: ", length(case_ids_joint), " samples")
  message("Joint ", fluid, " ctrl list saved: ", length(ctrl_ids_joint), " samples")

  # ------------------------------------------------------------------
  # 2. Load joint matrix (proteins x samples, NPX values = log2 scale)
  # ------------------------------------------------------------------
  mat_file <- file.path(path_matrix, paste0(fluid, "_joint.txt"))
  tmp      <- read.table(mat_file, header = TRUE, check.names = FALSE,
                         row.names = 1, sep = "\t", nrow = 10, quote = "")
  classes  <- sapply(tmp, class)
  mat      <- read.table(mat_file, header = TRUE, check.names = FALSE,
                         row.names = 1, sep = "\t", quote = "", colClasses = classes)

  colnames(mat) <- gsub('"', '', trimws(colnames(mat)))  # strip embedded quotes and whitespace
  rownames(mat) <- gsub('"', '', trimws(rownames(mat)))
  # ------------------------------------------------------------------
  # 3. Filter lists to samples actually present in the joint matrix
  # ------------------------------------------------------------------
  mat_samples <- colnames(mat)
  case_ids    <- case_ids_joint[case_ids_joint %in% mat_samples]
  ctrl_ids    <- ctrl_ids_joint[ctrl_ids_joint %in% mat_samples]

  message("After filtering to matrix columns:")
  message("  T2D (case): ", length(case_ids), " | Control: ", length(ctrl_ids))

  # ------------------------------------------------------------------
  # 4. Extract case / ctrl sub-matrices
  # ------------------------------------------------------------------
  dataCase <- mat[, case_ids, drop = FALSE]
  dataCtrl <- mat[, ctrl_ids, drop = FALSE]
  data_all  <- cbind(dataCase, dataCtrl)

  N <- ncol(dataCase)  # number of T2D samples
  M <- ncol(dataCtrl)  # number of control samples

  # NPX data is already on log2 scale — do NOT apply log2() again.

  # ------------------------------------------------------------------
  # 5. Remove proteins with >50% missing values
  # # ------------------------------------------------------------------
  # frac_na  <- rowMeans(is.na(data_all))
  # keep_na  <- which(frac_na <= 0.5)
  # dataCase <- dataCase[keep_na, , drop = FALSE]
  # dataCtrl <- dataCtrl[keep_na, , drop = FALSE]
  # data_all <- data_all[keep_na, , drop = FALSE]
  # message("Proteins after NA filter (<=50% missing): ", nrow(data_all))
  # 
  # # ------------------------------------------------------------------
  # # 6. Low-variation filter (bottom prc by IQR)
  # # ------------------------------------------------------------------
  # variation  <- apply(data_all, 1, IQR, type = 5, na.rm = TRUE)
  # soglia_prc <- quantile(variation, prc)
  # keep_iqr   <- which(variation > soglia_prc)
  # dataCase   <- dataCase[keep_iqr, , drop = FALSE]
  # dataCtrl   <- dataCtrl[keep_iqr, , drop = FALSE]
  # data_all   <- data_all[keep_iqr, , drop = FALSE]
  # message("Proteins after IQR filter (bottom ", prc * 100, "% removed): ", nrow(data_all))

  # ------------------------------------------------------------------
  # 7. Fold-change filter
  # log2FC = mean(T2D NPX) - mean(Ctrl NPX)  (NPX already on log2 scale)
  # ------------------------------------------------------------------
  logFC   <- rowMeans(dataCase, na.rm = TRUE) - rowMeans(dataCtrl, na.rm = TRUE)
  keep_fc  <- which(abs(logFC) >= soglia_log2fc)
  dataCase <- dataCase[keep_fc, , drop = FALSE]
  dataCtrl <- dataCtrl[keep_fc, , drop = FALSE]
  data_all <- data_all[keep_fc, , drop = FALSE]
  logFC    <- logFC[keep_fc]
  message("Proteins after FC filter (|log2FC| >= ", soglia_log2fc, "): ", length(logFC))

  # ------------------------------------------------------------------
  # 8. Unpaired Welch t-test (T2D vs Control are different subjects)
  # ------------------------------------------------------------------
  pval <- apply(data_all, 1, function(x) {
    x_case <- x[seq_len(N)]
    x_ctrl <- x[(N + 1):(N + M)]
    x_case <- x_case[!is.na(x_case)]
    x_ctrl <- x_ctrl[!is.na(x_ctrl)]
    if (length(x_case) < 3 || length(x_ctrl) < 3) return(NA_real_)
    t.test(x_case, x_ctrl, paired = FALSE)$p.value
  })

  # ------------------------------------------------------------------
  # 9. FDR correction (Benjamini-Hochberg)
  # ------------------------------------------------------------------
  pval_adj <- p.adjust(pval, method = "fdr")

  # ------------------------------------------------------------------
  # 10. Volcano plot (all proteins passing IQR + FC filter)
  # ------------------------------------------------------------------
  dirOut <- file.path(dirRes, fluid)
  if (!dir.exists(dirOut)) dir.create(dirOut, recursive = TRUE)

  sig_flag <- !is.na(pval_adj) & pval_adj <= soglia_pval

  fdr05_flag <- !is.na(pval_adj) & pval_adj <= 0.05

  # Build data frame for ggplot volcanos
  vdf <- data.frame(
    protein  = names(logFC),
    logFC    = logFC,
    neg_log10_fdr = -log10(pval_adj),
    fdr05    = fdr05_flag,
    stringsAsFactors = FALSE
  )

  # Asse x simmetrico attorno allo 0
  x_lim <- max(abs(vdf$logFC), na.rm = TRUE) * 1.05

  # --- Volcano 1: only FDR <= 0.05 proteins labeled in bold ---
  p1 <- ggplot(vdf, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(aes(color = fdr05), size = 1.5) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    # geom_vline(xintercept = c(-soglia_log2fc_plot, soglia_log2fc_plot),
    #            linetype = "dashed", color = "blue") +
    geom_text_repel(
      data     = subset(vdf, fdr05),
      aes(label = protein),
      fontface = "bold", size = 3,
      max.overlaps = Inf,
      box.padding = 0.4, point.padding = 0.2
    ) +
    scale_x_continuous(limits = c(-x_lim, x_lim)) +
    labs(title = paste0("Volcano: T2D vs Control (", fluid, ")"),
         x = "log2FC (T2D vs Control)", y = "-log10(FDR)") +
    theme_bw()
  ggsave(file.path(dirOut, "volcano.pdf"), plot = p1, width = 8, height = 7)

  # --- Volcano 2: all proteins labeled (base R, etichette sovrapposte) ---
  pdf(file.path(dirOut, "volcano_all_labels.pdf"), width = 10, height = 9)
  plot(vdf$logFC, vdf$neg_log10_fdr,
       main = paste0("Volcano (all labels): T2D vs Control (", fluid, ")"),
       xlab = "log2FC (T2D vs Control)",
       ylab = "-log10(FDR)",
       pch  = 16, cex = 0.7,
       col  = ifelse(vdf$fdr05, "red", "grey60"),
       xlim = c(-x_lim, x_lim))
  abline(h = -log10(0.05), lty = 2, lwd = 2, col = "red")
  # geom_vline(xintercept = c(-soglia_log2fc_plot, soglia_log2fc_plot),
  #            linetype = "dashed", color = "blue")
  text(vdf$logFC, vdf$neg_log10_fdr,
       labels = vdf$protein,
       cex    = 0.4, pos = 3,
       font   = ifelse(vdf$fdr05, 2, 1))
  dev.off()

  # ------------------------------------------------------------------
  # 11. Filter to significant proteins
  # ------------------------------------------------------------------
  keep_p    <- which(sig_flag)
  data_filt <- data_all[keep_p, , drop = FALSE]
  logFC_sig <- logFC[keep_p]
  pval_sig  <- pval_adj[keep_p]
  message("Significant proteins (FDR <= ", soglia_pval, "): ", length(keep_p))

  # ------------------------------------------------------------------
  # 12. Export DEG results table
  # ------------------------------------------------------------------
  direction <- ifelse(logFC_sig > 0, "up_in_T2D", "down_in_T2D")
  results   <- data.frame(
    protein   = rownames(data_filt),
    FDR       = pval_sig,
    log2FC    = logFC_sig,
    direction = direction,
    stringsAsFactors = FALSE
  )
  results <- results[order(results$log2FC, decreasing = TRUE), ]

  write.table(results,
              file      = file.path(dirOut, "DEG.txt"),
              row.names = FALSE, sep = "\t", quote = FALSE)

  message("DEG table saved: ", file.path(dirOut, "DEG.txt"))

  # ------------------------------------------------------------------
  # 13. Heatmap of significant proteins
  # ------------------------------------------------------------------
  if (nrow(data_filt) >= 2) {
    # Impute residual NAs with row mean before clustering
    data_heatmap <- t(apply(data_filt, 1, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    }))

    annotation <- data.frame(Group = c(rep("T2D", N), rep("Control", M)))
    rownames(annotation) <- colnames(data_filt)

    annotation_colors <- list(Group = c("T2D" = "orange", "Control" = "steelblue"))

    pheatmap(data_heatmap,
             scale                    = "row",
             border_color             = NA,
             cluster_cols             = TRUE,
             cluster_rows             = TRUE,
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             clustering_method        = "average",
             annotation_col           = annotation,
             annotation_colors        = annotation_colors,
             color                    = colorRampPalette(c("blue", "black", "yellow"))(100),
             show_rownames            = FALSE,
             show_colnames            = FALSE,
             cutree_cols              = 2,
             cutree_rows              = 2,
             width                    = 7,
             height                   = 7,
             filename                 = file.path(dirOut, "heatmap.pdf"))
    message("Heatmap saved: ", file.path(dirOut, "heatmap.pdf"))
  } else {
    message("Not enough significant proteins for heatmap (need >= 2).")
  }

  message("Done: ", fluid, " | Results in: ", dirOut)

  rm(mat, dataCase, dataCtrl, data_all, data_filt)
}

message("\nAll analyses complete.")
