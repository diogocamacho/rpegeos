#' rPEGEOS wrapper
#'
#' This function is a wrapper for a start to finish pathway enrichment analysis using tf-idf.
#' For a given gene set of interest, the algorithm will compute the cosine similarity between this vector and the gene sets in a given collection of pathways (for example, KEGG pathways.)
#' Additionally, the algorithm will generate a set of random gene sets that will be used to determine the probability of a given enrichment being random.
#' The output of the wrapper is a tidy data frame with all the enrichment scores that are better than 0.
#'
#' @param gene_set Query gene set of interest.  IDs on gene set need to match colnames of tf-idf matrix.
#' @param tfidf_matrix tf-idf matrix for pathway set. Column names are gene ids (EntrezIDs or gene symbols). Computed with \code{\link{tfidf}}
#' @param pathway_names Pathway names for pathway set, correspond to the rows of the tf-idf matrix.  Can be a Nx1 or Nx2 object: if Nx1, use pathway names; if Nx2, column 1 is pathway name, column 2 is pathway set source (e.g., KEGG, reactome, GO)
#' @param tfidf_crossproduct Cross-product of the tf-idf matrix. Computed with \code{\link{crossprod_matrix}}.
#' @param num_random Number of random sets to be generated to calculate significance of enrichment.  Defaults to 10,000.
#' @return A randomized set of gene sets as a matrix RxN where R is the number of random sets and N is the number of possible genes.
enrich_geneset <- function(gene_set, pathway_names, tfidf_matrix, tfidf_crossproduct, num_random)
{
  # run checks
  if(missing(gene_set)) stop("Need gene set.")
  if(length(gene_set) == 0) stop("No genes in gene set.")
  if(missing(num_random)) num_random <- 1000

  # match genes
  genes_per_pathway <- pathway_genes(gene_set = gene_set, pathway_tfidf = tfidf_matrix)

  # clean up pathways to match to:
  # we will remove those pathways that have no genes mapped to it
  nix <- which(genes_per_pathway$number_genes == 0)
  if(length(nix) != 0) {
    genes_per_pathway <- genes_per_pathway[-nix, ]
    tfidf_matrix <- tfidf_matrix[-nix, ]
    pathway_names <- pathway_names[-nix, ]
    tfidf_crossproduct <- tfidf_crossproduct[-nix]
  }

  # generate query vector
  # query_vector <- matrix(0,ncol=ncol(sparse_tfidf),nrow=1)
  query_vector <- integer(length = ncol(tfidf_matrix))
  query_vector[which(colnames(tfidf_matrix) %in% gene_set)] <- 1

  # cosine similarity between gene set and pathway_sets
  query_similarities <- cosine_similarity(tfidf_matrix = tfidf_matrix,
                                            query_vector = query_vector,
                                            tfidf_crossprod_mat = tfidf_crossproduct)

  # random probabilities
  prandom <- random_probability(similarity_results = query_similarities,
                                gs_size = sum(query_vector),
                                num_sets = num_random,
                                target_tfidf = tfidf_matrix,
                                tfidf_crossprod_mat = tfidf_crossproduct)

  # results
  res <- data_frame(cosine_similarity = as.vector(query_similarities),
                    probability_random = prandom,
                    number_genes = genes_per_pathway$number_genes,
                    genes_matched = genes_per_pathway$gene_names)

  # add pathway names and pathway source
  res <- res %>% tibble::add_column(., pathway_name = pathway_names, .before = 1)

  # put it all together
  res <- enrichment_results(query_similarities = query_similarities,
                            random_probability = prandom,
                            number_genes = genes_per_pathway$number_genes,
                            gene_names = genes_per_pathway$gene_names,
                            pathway_names = pathway_names)

  return(res)

}
