exportNetwork <- function(){
  
  ########################################
  # input parameters
  
  power <- input_parameter$power
  TOMType <- input_parameter$TOMType
  
  module <- input_parameter$module
  
  filename_dataInput <- output_file$filename_dataInput
  
  filename_networkConstruction <- output_file$filename_networkConstruction
  
  filename_cytoscapeInput_edges <-   output_file$filename_cytoscapeInput_edge
  filename_cytoscapeInput_nodes <- output_file$filename_cytoscapeInput_nodes
  ########################################
  
  lnames = load(file = filename_dataInput)
  
  lnames = load(file = filename_networkConstruction)
  ########################################   
  
  TOM = TOMsimilarityFromExpr(datExpr, power = power, TOMType = TOMType); 
  
  #modules = module
  modules = unique(moduleColors)
  
  nodes = names(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modNodes = nodes[inModule];
  
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modNodes, modNodes)
  ######################################## 
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = filename_cytoscapeInput_edges,
                                 nodeFile = filename_cytoscapeInput_nodes,
                                 weighted = TRUE,
                                 threshold = 0.05,
                                 nodeNames = modNodes,
                                 nodeAttr = moduleColors[inModule])

  # Add isolated nodes (no edge above threshold) as self-loops so they
  # appear in Cytoscape when importing only the edge file.
  edges <- read.table(filename_cytoscapeInput_edges, sep = "\t", header = TRUE,
                      quote = "", check.names = FALSE)
  nodes_in_edges <- unique(c(as.character(edges[, 1]), as.character(edges[, 2])))
  isolated       <- setdiff(modNodes, nodes_in_edges)

  if (length(isolated) > 0) {
    self_loops <- data.frame(
      fromNode = isolated,
      toNode   = isolated,
      weight   = 0,
      direction = "undirected",
      fromAltName = isolated,
      toAltName   = isolated,
      stringsAsFactors = FALSE
    )
    # Match column names from the exported edge file
    names(self_loops) <- names(edges)[1:ncol(self_loops)]
    write.table(self_loops, filename_cytoscapeInput_edges,
                sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE,
                col.names = FALSE)
    cat("Added", length(isolated), "isolated nodes as self-loops in edge file.\n")
  }

}





