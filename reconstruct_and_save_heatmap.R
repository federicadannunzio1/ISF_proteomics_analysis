# ==========================================
# 1. CARICAMENTO LIBRERIE E DATI
# ==========================================
library(WGCNA)

# Definiamo i percorsi (modifica se necessario)
dataset <- "SERUM" #"IF"
# Inserisci qui i percorsi ai file generati dalla Parte 1
base_dir <- paste0("~/Documents/projects/ISF/ISF_fede/code/WGCNA_paola/code/project/ISF/dataset/", dataset, "/Results/txtFile/")
file_cor  <- paste0(base_dir, "moduleTraitCor.txt")
file_pval <- paste0(base_dir, "moduleTraitPvalue.txt")

# Carichiamo le matrici
# check.names = FALSE serve a non trasformare gli spazi in punti nei nomi delle colonne
modTraitCor    <- as.matrix(read.table(file_cor, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
modTraitPvalue <- as.matrix(read.table(file_pval, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))

# ==========================================
# 2. PREPARAZIONE DEL TESTO (Stelline + Valori)
# ==========================================

# Creiamo le stelline basandoci sui p-value della tabella (logica grafica)
# * p < 0.05, ** p < 0.01, *** p < 0.001
stars <- ifelse(modTraitPvalue <= 0.001, "***", 
                ifelse(modTraitPvalue <= 0.01, "**", 
                       ifelse(modTraitPvalue <= 0.05, "*", "")))

# Componiamo la matrice di testo che apparirà nelle celle
# Formato: Correlazione* (a capo) (P-value)
textMatrix <- paste(signif(modTraitCor, 2), stars, "\n(",
                    signif(modTraitPvalue, 1), ")", sep = "")

# Fondamentale: paste trasforma tutto in un vettore, dobbiamo ridargli la forma di matrice
dim(textMatrix) <- dim(modTraitCor)

# ==========================================
# 3. GENERAZIONE HEATMAP
# ==========================================

output_pdf <- paste0("~/Documents/projects/ISF/ISF_fede/code/WGCNA_paola/code/project/ISF/dataset/", dataset, "/Results/figure/Heatmap_traits.pdf")

# Apriamo il dispositivo PDF (dimensioni 10x8 pollici, puoi variarle)
pdf(output_pdf, width = 12, height = 9)

# Fix Margini: c(bottom, left, top, right)
# Aumentiamo molto il primo valore (12 o 15) se i nomi dei tratti in basso sono lunghi
par(mar = c(15, 10, 3, 3))

# Creazione della heatmap
labeledHeatmap(
  Matrix = modTraitCor,
  xLabels = colnames(modTraitCor),   # Nomi dei tratti (Asse X)
  yLabels = rownames(modTraitCor),   # Nomi dei moduli (Asse Y)
  ySymbols = rownames(modTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),         # Scala colore: Blu (negativo), Bianco (0), Rosso (positivo)
  textMatrix = textMatrix,           # Testo con correlazione, stelle e p-val
  setStdMargins = FALSE,             # Usa i margini impostati con par(mar)
  cex.text = 0.6,                    # Dimensione del testo dentro le celle
  cex.lab.x = 0.8,                   # Dimensione etichette asse X
  cex.lab.y = 0.8,                   # Dimensione etichette asse Y
  zlim = c(-1, 1),                   # Range della scala colore
  main = "Module-trait relationships"
)

# Chiudiamo il file PDF
dev.off()

# Messaggio di conferma in console
cat("\n--- Analisi completata ---\n")
cat("Heatmap salvata in:", output_pdf, "\n")
