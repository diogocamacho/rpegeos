#' Significance of similarity
#'
#' \code{sig_sim} returns the probability of a given similarity score being better than random.
#'
#' @param similarity_results Sparse tf-idf matrix (computed from \code{\link{cosine_similarity}})
#' @param random_results A vector to compare to the tf-idf matrix. Length of vector is equal to number of columns of tf-idf matrix
#' @return Returns a vector with a probability of a given score being random
sig_sim <- function(similarity_results, random_results) {
  
  res_stats <- sapply(seq_along(similarity_results), 
                      function(x, y, i) length(which(y[i] > x[i, ]))/length(x[i,]), y = as.matrix(similarity_results) ,x = random_results)
  
  res_stats <- 1 - res_stats
  return(res_stats)
  
}