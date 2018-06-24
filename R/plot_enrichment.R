#' Plotting tf-idf enrichments
#'
#' Plots the tf-idf enrichments, ranked by similarity score.
#'
#' @param enrichment_results Results data frame
#' @param random_thr Threshold for probability of enrichment being random.  Defaults to 0.05.
#' @param max_pathways Number of top pathways to display.  Defaults to 25.
plot_enrichment <- function(results_df, random_thr, max_pathways)
{
  
  if(missing(max_pathways)) max_pathways <- 25
  if(missing(random_thr)) random_thr <- 0.05
  
  if (length(unique(results_df$pathway_source)) > 1) {
    results_df %>% 
      dplyr::filter(., probability_random < random_thr) %>% 
      dplyr::arrange(., desc(cosine_similarity)) %>% 
      dplyr::slice(1:max_pathways) %>% 
      ggplot(aes(x = fct_reorder(pathway_name, cosine_similarity), y = cosine_similarity)) + 
      geom_point(aes(color = pathway_source), size = 4) + 
      coord_flip() + 
      labs(y = "Similarity",x = "Pathway") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 20, face = "bold", color = "black"),
            axis.text = element_text(size = 12, color = "black"))
  } else {
    results_df %>% 
      dplyr::filter(., probability_random < random_thr) %>% 
      dplyr::arrange(., desc(cosine_similarity)) %>% 
      dplyr::slice(1:max_pathways) %>% 
      ggplot(aes(x = fct_reorder(pathway_name, cosine_similarity), y = cosine_similarity)) + 
      geom_point(size = 4) + 
      coord_flip() + 
      labs(y = "Similarity",x = "Pathway") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 20, face = "bold", color = "black"),
            axis.text = element_text(size = 12, color = "black"),
            legend.position = "none")
  }

}