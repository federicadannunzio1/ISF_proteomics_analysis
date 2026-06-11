rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("~/Documents/projects/ISF/ISF_fede/code/WGCNA/code/")

library(WGCNA)

# Source all functions directly (getSource() sources inside a function scope
# which doesn't expose symbols to the global env reliably)
source("src/script/input/config.R")
source("src/script/input/inputFiles.R")
source("src/script/output/outputDir.R")
source("src/script/output/outputFiles.R")
source("src/script/lib/dataInput.R")
source("src/script/lib/networkConstruction.R")
source("src/script/lib/relateModstoExt.R")
source("src/script/lib/visualization.R")
source("src/script/lib/exportNetwork.R")

# Parameters from parameters_wgcna.xlsx
# IF:   power=5, TOMType="unsigned"
# SERUM: power=4, TOMType="unsigned"

datasets <- list(
  IF    = list(dataset = "IF",    power = 5, TOMType = "unsigned",
               corType = "bicor", networkType = "unsigned",
               minModuleSize = 10, abline_h = 20,
               trait_name = "T2D", module = "turquoise"),
  SERUM = list(dataset = "SERUM", power = 4, TOMType = "unsigned",
               corType = "bicor", networkType = "unsigned",
               minModuleSize = 6,  abline_h = 11.5,
               trait_name = "T2D", module = "turquoise")
)

for (ds_name in names(datasets)) {

  cat("\n========== exportNetwork:", ds_name, "==========\n")

  params <- datasets[[ds_name]]
  path   <- paste0("project/ISF/dataset/", params$dataset)

  input_parameter <- list(
    path          = path,
    project       = "ISF",
    dataset       = params$dataset,
    filename_data   = paste0(path, "/matrix/matrix.txt"),
    filename_traits = paste0(path, "/matrix/Traits.txt"),
    abline_h      = params$abline_h,
    corType       = params$corType,
    power         = params$power,
    networkType   = params$networkType,
    TOMType       = params$TOMType,
    minModuleSize = params$minModuleSize,
    trait_name    = params$trait_name,
    module        = params$module
  )

  input_file  <- inputFiles()
  output_file <- outputFiles()

  exportNetwork()

  cat("Done:", ds_name, "\n")
}

cat("\nAll exports complete.\n")
