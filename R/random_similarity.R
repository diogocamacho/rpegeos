#' Random similarities
#'
#' \code{random_similarity} computes the cosine similarity between the set of randomized gene sets defined in \link{\code{random_sets}} and the pathway tf-idf matrix. It calls the \code{\link{cosine_similarity}} function.
#'
#' As the output of function we get a matrix that is PxR, where P is the number of pathways and R is the number of randomized sets.  This implies that each *column* in the matrix corresponds to the cosine similarity between that random set and each one of the pathways.
#'
#' @param random_set Matrix of randomized gene sets as defined on \code{\link{random_sets}}
#' @param tfidf_matrix Pathway tf-idf matrix.
#' @param tfidf_crossprod_mat Cross-product matrix for pathway tf-idf, computed with \code{\link{crossprod_matrix}}
#' @return A similarity matrix that is PxR, where P is the number of pathways and R is the number of random sets.
random_similarity <- function(random_set, tfidf_matrix, tfidf_crossprod_mat) {

  # rnd_sim <- apply(random_set,
  #                  1,
  #                  function(x) cosine_similarity(query_vector = x,
  #                                                tfidf_matrix = tfidf_matrix,
  #                                                tfidf_crossprod_mat = tfidf_crossprod_mat))

  rnd_sim <- rapply(object = random_set, f = cosine_similarity, tfidf_matrix = tfidf_matrix, tfidf_crossprod_mat = tfidf_crossprod_mat, how = "list")
  rnd_sim <- do.call(what = cbind, args = rnd_sim)


  # rnd_sim <- lapply(rnd_sim, FUN = as.matrix)
  # rnd_sim <- do.call(what = cbind, args = rnd_sim)

  return(rnd_sim)
}
