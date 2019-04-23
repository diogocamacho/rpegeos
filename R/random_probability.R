#' Random probability
#'
#' This function computes the cosine similarity between the set of randomized gene sets and the pathway tf-idf matrix. It calls the cosine_similarity function and returns the probability of a given cosine similarity being random.
#'
#' @param similarity_results Vector of cosine similarities
#' @param gene_set Query gene set
#' @param fthr Fold-change threshold for differential expression
#' @param pthr P-value threshold for differential expression
#' @param num_sets Number of random sets to generate
#' @param target_tfidf Pathway tf-idf matrix.
#' @param tfidf_crossprod_mat Cross-product matrix for pathway tf-idf.
#' @return A vector of random probabilities for each pathway given the gene set size.
random_probability <- function(similarity_results, gene_set, fthr, pthr, num_sets, target_tfidf, tfidf_crossprod_mat) {

  y1 <- matrix(0, ncol = ncol(target_tfidf), nrow = num_sets)
  y2 <- which(colnames(target_tfidf) %in% gene_set[, 1])
  x2 <- length(with(gene_set, which(abs(gene_set[, 2]) > fthr & gene_set[, 3] < pthr)))

  # this generates the random set from the query data (!) before computing similarity
  for (i in seq(1, nrow(y1))) {
    a <- sample(y2, size = x2, replace = FALSE)
    y1[i, a] <- 1
  }

  # now compute cosine similarity
  a1 <- apply(y1, 1, crossprod)
  a2 <- sapply(a1, function(y) sqrt(cpm * y))
  a3 <- apply(y1, 1, function(y) tfidf_matrix %*% y)
  cs <- a3 / a2 # <-- cosine similarity

  # compute probability of being random
  ps <- matrix(0, nrow = nrow(a3), ncol = ncol(a3))
  for (i in seq(1, nrow(ps))) {
    a <- which(cs[i, ] > similarity_results[i])
    ps[i, a] <- 1
  }
  ps <- rowSums(ps) / ncol(ps)

  return(ps)
}
