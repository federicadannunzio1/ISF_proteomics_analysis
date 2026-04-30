config <- function(){
  
  #################################
  #project <- "TutorialWGCNA"
  #dataset <- "FemaleLiver3600"
  
  project <- "ISF"
  dataset <-  "SERUM" #"IF"
  
  path <- paste0("project/",project,"/dataset/",dataset)
  #################################
  # input files
  
  filename_data <- paste0(path,"/matrix/matrix.txt")
  
  filename_traits <- paste0(path,"/matrix/Traits.txt")
  #################################
  # input parameters
  
  abline_h <- 11.5
  
  corType <- "bicor" 
  power <- 4
  networkType <-  "unsigned"
  TOMType <- "unsigned"
  minModuleSize <- 6
  
  trait_name <- "T2D"
  
  module <- "turquoise"
  #################################
  
  input_parameter <- list(path = path,
                          project = project,
                          dataset = dataset,
                          filename_data = filename_data,
                          filename_traits = filename_traits,
                          abline_h = abline_h,
                          corType = corType,
                          power = power,
                          networkType = networkType,
                          TOMType = TOMType,
                          minModuleSize = minModuleSize,
                          trait_name = trait_name,
                          module = module)
  
  return(input_parameter)
  
}