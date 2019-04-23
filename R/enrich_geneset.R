#' rPEGEOS wrapper
#'
#' This function is a wrapper for a start to finish pathway enrichment analysis using tf-idf.
#' For a given gene set of interest, the algorithm will compute the cosine similarity between this vector and the gene sets in a given collection of pathways (for example, KEGG pathways.)
#' Additionally, the algorithm will generate a set of random gene sets that will be used to determine the probability of a given enrichment being random.
#' The output of the wrapper is a tidy data frame with all the enrichment scores that are better than 0.
#'
#' @param gene_set Query gene set of interest. Should be a Nx3 matrix, where column 1 is Entrez ID, column 2 is fold-change, and column 3 is p-value (recommend using FDR corrected p-values from a package like `limma` or `DESeq2`)
#' @return A tibble with columns geneset, number_genes, cosine_similarity, and probability_random.
enrich_geneset <- function(gene_set)
{

  message("Checks and balances...")
  if(missing(gene_set)) stop("Need gene set.")
  if(nrow(gene_set) == 0) stop("No genes in gene set.")
  if(ncol(gene_set) != 3) stop("Incorrect format for gene set. Gene set is an Nx3 matrix where column 1 is EntrezID, column 2 is fold change, and column 3 is fold change p-value.")

  message("--- Pathway enrichments with rPEGEOS ---")
  message("")

  gs <- geneset_selection()
  ui <- user_inputs()

  tmp <- rpegeos::pathway_sets[[gs]]

  # match genes
  message("Counting gene representation...")
  genes_per_pathway <- pathway_genes(gene_set = gene_set,
                                     pathway_tfidf = tmp$tfidf,
                                     fthr = ui[[1]],
                                     pthr = ui[[2]])

  # clean up pathways to match to:
  # we will remove those pathways that have no genes mapped to it
  message("Cleaning up gene sets...")
  nix <- which(genes_per_pathway$num_diff == 0)
  if(length(nix) != 0) {
    genes_per_pathway <- genes_per_pathway[-nix, ]
    tfidf_matrix <- tmp$tfidf[-nix, ]
    pathway_names <- tmp$geneset_name[-nix]
    cpm <- tmp$cpm[-nix]
  }

  # generate query vector
  # query_vector <- matrix(0,ncol=ncol(sparse_tfidf),nrow=1)
  message("Generating query vector...")
  query_vector <- integer(length = ncol(tfidf_matrix))
  query_vector[which(colnames(tfidf_matrix) %in% gene_set[, 1])] <- 1

  # cosine similarity between gene set and pathway_sets
  message("Computing cosine similarities...")
  query_similarities <- cosine_similarity(tfidf_matrix = tfidf_matrix,
                                            query_vector = query_vector,
                                            tfidf_crossprod_mat = cpm)

  # random probabilities
  message("Running random sets...")
  prandom <- random_probability(similarity_results = query_similarities,
                                gene_set = gene_set,
                                fthr = ui[[1]],
                                pthr = ui[[2]],
                                num_sets = ui[[3]],
                                target_tfidf = tfidf_matrix,
                                tfidf_crossprod = cpm)

  # results
  message("Compiling results...")
  res <- enrichment_results(query_similarities = query_similarities,
                            random_probability = prandom,
                            gene_counts = genes_per_pathway,
                            pathway_names = pathway_names)

  res <- res %>%
    dplyr::mutate(., probability_random = replace(probability_random, probability_random == 0, 1/ui[[3]])) %>%
    dplyr::mutate(., enrichment_score = cosine_similarity - log10(probability_random) + 1) %>%
    dplyr::arrange(., desc(enrichment_score))

  # plotting results
  # first, volcano plot of data
  # p1 <- plot_volcano(gene_set = gene_set,
  #                    fthr = ui[[1]],
  #                    pthr = ui[[2]])

  # now enrichment results
  # p2 <- plot_enrichment(results_df = res)
  # p2

  return(res)
  message("")
  message("----------")
  message("Summary: ")
  message(paste("Number of pathways with probability random < 0.001: ",
                res %>% dplyr::filter(., probability_random < 0.001) %>% nrow))
  message(paste("Number of pathways with probability random < 0.01: ",
                res %>% dplyr::filter(., probability_random < 0.01) %>% nrow))
  message(paste("Top 10 enriched pathways: "))
  res %>% dplyr::slice(1:10) %>% dplyr::select(., geneset, enrichment_score)
  message("----------")
  message("")
  message("** DONE **")
}
