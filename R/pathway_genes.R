#' Genes in pathway
#'
#' Returns the number of genes in the query set that are present in a given pathway.
#'
#' @param gene_set Query gene set.  IDs of genes match colnames of tf-idf matrix.
#' @param pathway_tfidf tf-idf matrix computed with \code{\link{tfidf}}
#' @return Returns vector with number of genes matched to any given pathway.
pathway_genes <- function(gene_set, pathway_tfidf) {

  x <- which(colnames(pathway_tfidf) %in% gene_set)
  y <- pathway_tfidf[, x]
  y[y != 0] <- 1
  y <- rowSums(y)

  return(y)

}
