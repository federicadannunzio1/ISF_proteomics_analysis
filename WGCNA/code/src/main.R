rm(list = ls())

options(stringsAsFactors = FALSE)

setwd("~/Documents/projects/ISF/ISF_fede/code/WGCNA/code/")
######################################
library(WGCNA)

source("src/script/getSource.R")
######################################
getSource()
input_parameter <- config()
input_file <- inputFiles()
output_file <- outputFiles()
######################################

dataInput()

networkConstruction()

relateModstoExt()

visualization()

exportNetwork()
