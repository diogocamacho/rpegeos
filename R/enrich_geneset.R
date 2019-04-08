#' rPEGEOS wrapper
#'
#' This function is a wrapper for a start to finish pathway enrichment analysis using tf-idf.
#' For a given gene set of interest, the algorithm will compute the cosine similarity between this vector and the gene sets in a given collection of pathways (for example, KEGG pathways.)
#' Additionally, the algorithm will generate a set of random gene sets that will be used to determine the probability of a given enrichment being random.
#' The output of the wrapper is a tidy data frame with all the enrichment scores that are better than 0.
#'
#' @param gene_set Query gene set of interest, provided as Entrez gene IDs.
#' @return A randomized set of gene sets as a matrix RxN where R is the number of random sets and N is the number of possible genes.
enrich_geneset <- function(gene_set)
{

  message("Checks and balances...")
  if(missing(gene_set)) stop("Need gene set.")
  if(length(gene_set) == 0) stop("No genes in gene set.")


  message("--- Pathway enrichments with rPEGEOS ---")
  message("")

  gs <- geneset_selection()
  message("")
  nrandom <- user_inputs()
  message("")

  tmp <- rpegeos::pathway_sets[[gs]]

  # match genes
  genes_per_pathway <- pathway_genes(gene_set = gene_set, pathway_tfidf = tmp$tfidf)

  # clean up pathways to match to:
  # we will remove those pathways that have no genes mapped to it
  nix <- which(genes_per_pathway == 0)
  if(length(nix) != 0) {
    genes_per_pathway <- genes_per_pathway[-nix]
    tfidf_matrix <- tmp$tfidf[-nix, ]
    pathway_names <- tmp$geneset_name[-nix]
    cpm <- tmp$cpm[-nix]
  }

  # generate query vector
  # query_vector <- matrix(0,ncol=ncol(sparse_tfidf),nrow=1)
  query_vector <- integer(length = ncol(tfidf_matrix))
  query_vector[which(colnames(tfidf_matrix) %in% gene_set)] <- 1

  # cosine similarity between gene set and pathway_sets
  query_similarities <- cosine_similarity(tfidf_matrix = tfidf_matrix,
                                            query_vector = query_vector,
                                            tfidf_crossprod_mat = cpm)

  # random probabilities
  prandom <- random_probability(similarity_results = query_similarities,
                                gs_size = sum(query_vector),
                                num_sets = nrandom,
                                target_tfidf = tfidf_matrix,
                                tfidf_crossprod_mat = cpm)

  # results
  res <- enrichment_results(query_similarities = query_similarities,
                            random_probability = prandom,
                            number_genes = genes_per_pathway,
                            pathway_names = pathway_names)

  return(res)

}
