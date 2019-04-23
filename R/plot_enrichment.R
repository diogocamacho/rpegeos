#' Plotting tf-idf enrichments
#'
#' Plots the tf-idf enrichments, ranked by similarity score.
#'
#' @param results_df Results data frame
#' @param max_pathways Number of top pathways to display.  Defaults to 20.
plot_enrichment <- function(results_df) {


  if (nrow(results_df) > 10) {
    results_df %>%
      dplyr::arrange(., desc(enrichment_score)) %>%
      dplyr::slice(1:10) %>%
      ggplot() +
      geom_point(aes(x = enrichment_score, y = forcats::fct_reorder(geneset, enrichment_score), size = number_genes)) +
      labs(x = "Enrichment score", y = NULL) +
      theme_bw() +
      theme(axis.text = element_text(size = 18, color = "black"),
            axis.title = element_text(size = 24, color = "black"),
            legend.position = "none",
            panel.grid = element_blank())
  } else {
    results_df %>%
      dplyr::arrange(., desc(enrichment_score)) %>%
      ggplot() +
      geom_point(aes(x = enrichment_score, y = forcats::fct_reorder(geneset, enrichment_score), size = number_genes)) +
      labs(x = "Enrichment score", y = NULL) +
      theme_bw() +
      theme(axis.text = element_text(size = 18, color = "black"),
            axis.title = element_text(size = 24, color = "black"),
            legend.position = "none",
            panel.grid = element_blank())
  }
}
