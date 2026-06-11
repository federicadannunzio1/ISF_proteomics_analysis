rm(list = ls())

library(tidyverse)
library(xlsx)

sample_map_raw <- read.xlsx("~/Documents/collaborations/ISF/ISF_fede/data/samples_info.xlsx", 
                            sheetIndex = 1)

colnames(sample_map_raw) <- c("ctrl_serum", "ctrl_if", "subject", "t2b_serum", "t2b_if", "subject")

sample_map_raw <- sample_map_raw %>% 
  select(c("ctrl_serum","ctrl_if","t2b_serum", "t2b_if"))

sample_map_raw <- sample_map_raw[-1,]
####################################

IF_ctrl_list <- sample_map_raw %>%
  select(ctrl_if) %>%
  filter(!is.na(ctrl_if)) %>%
  rename(SampleID = ctrl_if) %>%
  mutate(Group = "IF_ctrl") %>%
  select(SampleID)

serum_ctrl_list <- sample_map_raw %>%
  select(ctrl_serum) %>%
  filter(!is.na(ctrl_serum)) %>%
  rename(SampleID = ctrl_serum) %>%
  mutate(Group = "serum_ctrl") %>%
  select(SampleID)

IF_t2d_list <- sample_map_raw %>%
  select(t2b_if) %>%
  filter(!is.na(t2b_if)) %>%
  rename(SampleID = t2b_if) %>%
  mutate(Group = "IF_t2b") %>%
  select(SampleID)

serum_t2d_list <- sample_map_raw %>%
  select(t2b_serum) %>%
  filter(!is.na(t2b_serum)) %>%
  rename(SampleID = t2b_serum) %>%
  mutate(Group = "serum_t2b") %>%
  select(SampleID)
#################################
# Keeping in the list just the samples that are present in the matrix

all_lists <- list(IF_ctrl_list, IF_t2d_list, serum_ctrl_list, serum_t2d_list)

files <- list.files(
  "~/Documents/collaborations/ISF/ISF_fede/data/matrix",
  full.names = TRUE,
  pattern = "\\.txt$"
)

filtered_lists <- list()
dirRes <-  "~/Documents/collaborations/ISF/ISF_fede/data/list/"

for (f in files) {
  
  fname <- basename(f)
  
  mat <- read.table(
    f,
    sep = "\t",
    header = TRUE,
    quote = "",
    check.names = FALSE,
    row.names = 1
  )
  
  mat_samples <- trimws(colnames(mat))
  
  if (grepl("_IF_", fname)) {
    list_idx <- c(1, 2)   # IF ctrl, IF case
    list_names <- c("ctrl", "case")
  } else if (grepl("_serum_", fname)) {
    list_idx <- c(3, 4)   # Serum ctrl, Serum case
    list_names <- c("ctrl", "case")
  } else {
    stop("Tipo matrice non riconosciuto: ", fname)
  }

  for (j in seq_along(list_idx)) {
    
    df <- all_lists[[list_idx[j]]]
    df$SampleID <- trimws(df$SampleID)
    
    df_filt <- df[df$SampleID %in% mat_samples, , drop = FALSE]
    
    filtered_lists[[paste0(fname, "_", list_names[j])]] <- df_filt
    
    message(
      "Done: ", fname,
      " | ", list_names[j],
      " | SampleID trovati: ", nrow(df_filt)
    )
  }
  

  corrected_fname <- strsplit(fname, "\\.")[[1]][1]
  corrected_fname <- gsub("matrix", "list", corrected_fname)
  
  for (j in seq_along(list_names)) {
    
    list_key <- paste0(fname, "_", list_names[j])
    df_to_save <- filtered_lists[[list_key]]
    colnames(df_to_save) <- NULL
    
    write.table(
      df_to_save,
      file = paste0(dirRes, corrected_fname, "_", list_names[j], ".txt"),
      sep = "\t",
      quote = FALSE,
      row.names = F
    )
  }
  
}
