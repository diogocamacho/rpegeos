#' Results data frame
#'
#' This function returns a tidy data frame with all the enrichment results.
#'
#' @param query_similarities Similarity score for query vector
#' @param random_probability Probability of a given enrichment being random
#' @param number_genes Number of genes as described in \code{\link{pathway_genes}}
#' @param gene_names Gene names as described in \code{\link{pathway_genes}}
#' @param pathway_names Pathway names for pathway set, correspond to the rows of the tf-idf matrix.  Can be a Nx1 or Nx2 object: if Nx1, use pathway names; if Nx2, column 1 is pathway name, column 2 is pathway set source (e.g., KEGG, reactome, GO)
#' @return Returns the enrichment results for each pathway given the query set.
enrichment_results <- function(query_similarities, random_probability, number_genes, gene_names, pathway_names) {

  res <- data_frame(cosine_similarity = as.vector(query_similarities),
                    probability_random = random_probability,
                    number_genes = number_genes,
                    genes_matched = gene_names)

  # add pathway names and pathway source
  res <- res %>%
    tibble::add_column(., pathway_name = pathway_names, .before = 1)

  return(res)
}
