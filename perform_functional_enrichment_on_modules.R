rm(list = ls())
library(enrichR)
library(ggplot2)
library(xlsx)

# ---------------------------------------------------------------
# PARAMETRI GLOBALI
# ---------------------------------------------------------------
pval_thr    <- 0.05
# Qui puoi aggiungere tutti i DB che desideri (es. "KEGG_2021_Human", "WikiPathway_2021_Human")

dbs <- c("GO_Biological_Process_2025", "KEGG_2026", "GO_Molecular_Function_2025", "DisGeNET")
all_modules <- c("turquoise", "blue", "brown", "yellow", "green") # remember to add module green if I'm considering isf

# Cartella di output
dataset <-"IF" #"SERUM"
res_path <- paste0("~/Documents/projects/ISF/ISF_fede/code/WGCNA_paola/code/project/ISF/dataset/", dataset, "/Results/functional_enrichment")

if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
} else { 
  print("the dir already exists")
}

# Connessione a Enrichr
setEnrichrSite("Enrichr")

# ---------------------------------------------------------------
# CARICAMENTO DATI
# ---------------------------------------------------------------
input_path <- paste0("~/Documents/projects/ISF/ISF_fede/code/WGCNA_paola/code/project/ISF/dataset/", dataset, "/Results/txtFile")
gene_file <- file.path(input_path, "geneInfo.txt")
df <- read.delim(gene_file, header = TRUE, sep = "\t")

# ---------------------------------------------------------------
# LOOP ESTERNO: Moduli
# ---------------------------------------------------------------
for (mod in all_modules) {
  
  cat("\nAnalisi Modulo:", mod)
  
  # Selezione geni
  genes_to_test <- df$geneSymbol[df$moduleColor == mod]
  if(is.null(genes_to_test)) genes_to_test <- df[df$moduleColor == mod, 1] 
  
  genes_to_test <- na.omit(genes_to_test)
  
  if (length(genes_to_test) < 3) {
    cat("  -> Troppi pochi geni - Salto il modulo.\n")
    next
  }
  
  # Richiamo Enrichr per tutti i DB selezionati
  enriched <- tryCatch(
    enrichr(genes_to_test, dbs),
    error = function(e) return(NULL)
  )
  
  if (is.null(enriched)) {
    cat("  -> Errore di connessione a Enrichr - Salto il modulo.\n")
    next
  }
  
  # ---------------------------------------------------------------
  # LOOP INTERNO: Database (Novità: cicla sui DB richiesti)
  # ---------------------------------------------------------------
  for (current_db in dbs) {
    
    cat("\n  -> Database:", current_db)
    
    res_table <- enriched[[current_db]]
    
    # Funzione di filtraggio originale
    process_res <- function(res_table) {
      if (is.null(res_table) || !is.data.frame(res_table) || nrow(res_table) == 0) {
        return(NULL)
      }
      filtered <- res_table[res_table$Adjusted.P.value < pval_thr, ]
      if (nrow(filtered) == 0) return(NULL)
      return(filtered)
    }
    
    res_filtered <- process_res(res_table)
    
    if (is.null(res_filtered)) {
      cat("  -> Nessun termine significativo trovato per", current_db)
      next
    }
    
    # 1. Salvataggio Excel con nome DINAMICO del Database
    xlsx_file <- file.path(res_path, paste0(mod, "_", current_db, "_enrichment.xlsx"))
    write.xlsx(res_filtered, xlsx_file, row.names = FALSE)
    
    # 2. Preparazione Top 10 per Plot
    res_top10 <- head(res_filtered[order(res_filtered$Adjusted.P.value), ], 10)
    res_top10$Count <- as.numeric(sub("/.*", "", res_top10$Overlap))
    
    # 3. Plot con Titolo dinamico
    p <- ggplot(res_top10, aes(x = Combined.Score, y = reorder(Term, Combined.Score))) +
      geom_point(aes(size = Count, color = Adjusted.P.value)) +
      scale_color_gradient(low = "red", high = "blue") +
      scale_size_continuous(range = c(5, 12)) +
      theme_classic() +
      labs(
        title = paste("Top 10 terms -", current_db,"-", mod),
        subtitle = paste0("Module: ", mod, " | Adj. p-value < ", pval_thr),
        x = "Combined Score", y = ""
      )
    
    # 4. Salvataggio PDF con nome DINAMICO del Database
    pdf_file <- file.path(res_path, paste0(mod, "_", current_db, "_plot.pdf"))
    ggsave(filename = pdf_file, plot = p, width = 12, height = 7)
  }
}

cat("\n\nProcesso completato con successo.\n")