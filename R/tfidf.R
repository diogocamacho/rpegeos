#' TF-IDF matrix
#'
#' Returns a tf-idf sparse matrix given a document-term matrix.
#' This matrix is a DxT matrix, where D is the set of documents in a corpus and T are the unique terms in that corpus.
#' The matrix reflects a count of occurences of a given term in a given document.
#'
#'
#' @param data_matrix Sparse document-term matrix (generated with Matrix package)
#' @return A sparse matrix with the computed tf-idf
rpegeos_tfidf <- function(data_matrix)
{
  # tidy data ----
  dtm1 <- tidytext::tidy(data_matrix)

  # count words in dtm ----
  words1 <- dtm1 %>%
    dplyr::count(row, column) %>%
    dplyr::ungroup()

  # count all terms ----
  terms1 <- words1 %>%
    dplyr::group_by(row) %>%
    dplyr::summarize(total = sum(n))

  # calculate tf-idf ----
  tfidf1 <- words1 %>%
    tidytext::bind_tf_idf(column, row, n)

  # build matrix ----
  # output matrix is tfidf1 in same order as data_matrix

  m1 <- matrix(0, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  for (i in seq(1, nrow(data_matrix))) {
    a1 <- which(tfidf1$row == i)
    a2 <- which(colnames(data_matrix) %in% tfidf1$column[a1])
    m1[i, a2] <- tfidf1$tf_idf[a1]
  }
  colnames(m1) <- colnames(data_matrix)

  # # TIDY DATA FIRST
  # tidy_dtm <- tidytext::tidy(data_matrix)
  #
  # # document-term matrix processing
  # doc_words <- tidy_dtm %>%
  #   dplyr::count(row, column) %>%
  #   ungroup()
  #
  # total_terms <- doc_words %>%
  #   dplyr::group_by(row) %>%
  #   dplyr::summarize(total = sum(n))
  #
  # tfidf <- doc_words %>%
  #   tidytext::bind_tf_idf(column, row, n)
  #
  # tfidf <- tidytext::cast_sparse(data = tfidf, row, column, value = tf_idf)

  # return(tfidf)
  return(m1)
}
