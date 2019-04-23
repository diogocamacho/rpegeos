#' Results data frame
#'
#' This function returns a tidy data frame with all the enrichment results.
#'
#' @param query_similarities Similarity score for query vector
#' @param random_probability Probability of a given enrichment being random
#' @param gene_counts Number of genes as described in \code{\link{pathway_genes}}
#' @param pathway_names Pathway names for pathway set, correspond to the rows of the tf-idf matrix.
#' @return Returns the enrichment results for each pathway given the query set.
enrichment_results <- function(query_similarities, random_probability, gene_counts, gene_names, pathway_names) {

  res <- tibble::tibble(cosine_similarity = as.vector(query_similarities),
                        probability_random = random_probability,
                        number_genes = gene_counts[, 3])

  # add pathway names and pathway source
  res <- res %>%
    tibble::add_column(., geneset = pathway_names, .before = 1)

  return(res)
}
