rm(list = ls())

library(ggvenn)

if_res <- read.table("~/Documents/projects/ISF/ISF_fede/code/WGCNA_paola/code/project/ISF/dataset/IF/Results/txtFile/geneInfo.txt", sep = "\t", header = T, quote = "", check.names = F, 
                     row.names = 1)
rownames(if_res) <- gsub('"', '', rownames(if_res))

serum_res <- read.table("~/Documents/projects/ISF/ISF_fede/code/WGCNA_paola/code/project/ISF/dataset/SERUM/Results/txtFile/geneInfo.txt", sep = "\t", header = T, quote = "", check.names = F,
                        row.names = 1)
rownames(serum_res) <- gsub('"', '', rownames(serum_res))

if_turquoise_module <- if_res[which(if_res$moduleColor == "turquoise"),]
serum_brown_module <- serum_res[which(serum_res$moduleColor == "brown"),]

x <- intersect(rownames(if_turquoise_module), rownames(serum_brown_module))

list_for_venn <- list(rownames(if_turquoise_module), rownames(serum_brown_module))
names(list_for_venn) <- c("IF_turquoise_module", "SERUM_brown_module")

pdf("~/Documents/projects/ISF/ISF_fede/code/modules_analysis/venn.pdf")
ggvenn(list_for_venn,  set_name_color = c("turquoise", "brown"), stroke_color = c("turquoise", "brown"),
       fill_color = c("turquoise", "brown"), show_percentage = F)
dev.off()

max_ln <- max(length(x), 
              length(setdiff(rownames(if_turquoise_module), rownames(serum_brown_module))), 
              length(setdiff(rownames(serum_brown_module), rownames(if_turquoise_module))))

fill_vec <- function(v, len) {
  c(v, rep("", len - length(v)))
}

summary_table <- data.frame(
  common = fill_vec(x, max_ln),
  unique_in_if = fill_vec(setdiff(rownames(if_turquoise_module), rownames(serum_brown_module)), max_ln),
  unique_in_serum = fill_vec(setdiff(rownames(serum_brown_module), rownames(if_turquoise_module)), max_ln)
)

write.table(summary_table, "~/Documents/projects/ISF/ISF_fede/code/modules_analysis/summary_table_common_proteins.txt",
            sep = "\t",row.names = F, quote = F)

