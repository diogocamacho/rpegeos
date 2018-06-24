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
  # TIDY DATA FIRST
  tidy_dtm <- tidytext::tidy(data_matrix)

  # document-term matrix processing
  doc_words <- tidy_dtm %>%
    dplyr::count(row, column) %>%
    ungroup()

  total_terms <- doc_words %>%
    dplyr::group_by(row) %>%
    dplyr::summarize(total = sum(n))

  tfidf <- doc_words %>%
    tidytext::bind_tf_idf(column, row, n)

  tfidf <- tidytext::cast_sparse(data = tfidf, row, column, value = tf_idf)

  return(tfidf)
}
