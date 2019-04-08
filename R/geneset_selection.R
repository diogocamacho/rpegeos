#' Select gene set source.
#'
#' Function with no parameters, allows user to choose which gene set to use for enrichment.
geneset_selection <- function() {
  message("Available gene sets:")
  message("[1] KEGG")
  message("[2] PharmGKB")
  message("[3] Small Molecule Pathway DB")
  message("[4] WikiPathways")
  message("[5] SignaLink")
  message("[6] HumanCyc")
  message("[7] Biocarta")
  message("[8] Reactome")
  message("[9] INOH")
  message("[10] PID")
  message("[11] NetPath")
  message("[12] EHMN")
  message("[13] Gene Ontology: Biological Processes")
  message("[14] Gene Ontology: Cellular Component")
  message("[15] Gene Ontology: Molecular Function")
  message("[16] Broad Hallmark")
  message("[17] Run all")
  message("")
  selection <- readline(promp = "Selection: ")
  selection <- as.integer(selection)
  if (selection %in% seq(1, 17)) {
    return(selection)
  }
  else {
    message("Not an available data set. Stopping.")
  }
}
