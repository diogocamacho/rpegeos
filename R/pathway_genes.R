#' Genes in pathway
#'
#' Returns the number of genes in the query set that are present in a given pathway.
#'
#' @param gene_set Query gene set of interest. Should be a Nx3 matrix, where column 1 is Entrez ID, column 2 is fold-change, and column 3 is p-value (recommend using FDR corrected p-values from a package like `limma` or `DESeq2`)
#' @param pathway_tfidf tf-idf matrix computed with \code{\link{tfidf}}
#' @param fthr Fold change threshold
#' @param pthr p-value threshold
#' @return Returns vector with number of genes matched to any given pathway.
pathway_genes <- function(gene_set, pathway_tfidf, fthr, pthr) {


  up_genes <- gene_set[which(gene_set[, 2] > fthr & gene_set[, 3] < pthr), 1]
  down_genes <- gene_set[which(gene_set[, 2] < -fthr & gene_set[, 3] < pthr), 1]

  x_up <- which(colnames(pathway_tfidf) %in% up_genes)
  x_down <- which(colnames(pathway_tfidf) %in% down_genes)
  x_total <- which(colnames(pathway_tfidf) %in% gene_set[, 1]) # <-- possible number of genes in gene set
  path_count <- pathway_tfidf[, x_total]
  path_count[path_count != 0] <- 1
  path_count <- rowSums(path_count)

  if (length(x_up) > 1) {
    y_up <- pathway_tfidf[, x_up]
    y_up[y_up != 0] <- 1
    y_up <- rowSums(y_up)
  } else if (length(x_up) == 1) {
    y_up <- pathway_tfidf[, x_up]
    y_up[y_up != 0] <- 1
  } else {
    y_up <- integer(length = nrow(pathway_tfidf))
  }


  if (length(x_down) > 1) {
    y_down <- pathway_tfidf[, x_down]
    y_down[y_down != 0] <- 1
    y_down <- rowSums(y_down)
  } else if (length(x_down) == 1) {
    y_down <- pathway_tfidf[, x_down]
    y_down[y_down != 0] <- 1
  } else {
    y_down <- integer(length = nrow(pathway_tfidf))
  }


  gene_counts <- tibble::tibble(num_up = y_up,
                                num_down = y_down,
                                num_diff = y_up + y_down,
                                tot = path_count,
                                prop = (y_up + y_down) / path_count)

  return(gene_counts)

}
