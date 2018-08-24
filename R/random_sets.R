#' Generate random sets
#'
#' \code{random_sets} returns a set of randomized gene sets to assess significance of enriched pathways.
#'
#' @param number_sets How many random sets to generate.  Defaults to 10,000
#' @param universe_size Size of total possible number of genes.  User defined, should match number of columns in tf-idf matrix.
#' @param gene_set_size How many genes to sample.  Should be the same size as the query vector.
#' @return A randomized set of gene sets as a matrix RxN where R is the number of random sets and N is the number of possible genes.
random_sets <- function(number_sets, universe_size, gene_set_size) {

  if(missing(number_sets)) number_sets <- 10000
  if(missing(universe_size)) stop("Need to know universe size. Exiting.")
  if(missing(gene_set_size)) stop("Need to know size of gene set to sample. Exiting.")
  if(length(gene_set_size) == 0) stop("No genes in gene set. Exiting.")


  random_set <- vector(mode = "list", length = number_sets)
  for(i in seq(1, number_sets)) {
    s1 <- sample(universe_size, size = gene_set_size, replace = FALSE)
    v1 <- integer(universe_size)
    v1[s1] <- 1
    random_set[[i]] <- v1
  }
  # random_set <- Matrix::Matrix(data = do.call(what = rbind, args = random_set), sparse = TRUE)

  return(random_set)
}
