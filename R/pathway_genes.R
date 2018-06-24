#' Genes in pathway
#'
#' Returns the number and ids of the genes in the query set that are present in a given pathway.
#'
#' @param gene_set Query gene set.  IDs of genes match colnames of tf-idf matrix.
#' @param pathway_tfidf tf-idf matrix computed with \code{\link{rpegeos_tfidf}}
#' @return Returns Nx2 tibble with number of genes and gene names matched to any given pathway.
pathway_genes <- function(gene_set, pathway_tfidf) {

  num_genes <- integer(nrow(pathway_tfidf))
  gene_names <- rep(NA,nrow(pathway_tfidf))

  for(i in seq(1, nrow(pathway_tfidf))) {
    all_found <- intersect(names(which(pathway_tfidf[i,] != 0)),gene_set)
    num_genes[i] <- length(all_found)
    gene_names[i] <- paste(all_found,collapse = " | ")
  }

  gp <- data_frame(number_genes = num_genes,
                   gene_names = gene_names)

  return(gp)

}
