plot_volcano <- function(gene_set, fthr, pthr) {

  df <- tibble::tibble(gene = gene_set[, 1],
                       fold = gene_set[, 2],
                       pval = gene_set[, 3])

  df %>%
    tibble::add_column(., color = "black") %>%
    dplyr::mutate(., color = replace(color, list = which(fold > fthr & pval < pthr), values = "red")) %>%
    dplyr::mutate(., color = replace(color, list = which(fold < -fthr & pval < pthr), values = "red")) %>%
    ggplot() +
    geom_point(aes(x = fold, y = -log10(pval), color = color), alpha = 0.5) +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = -log10(pthr), lty = 2) +
    geom_vline(xintercept = -fthr, lty = 2) +
    geom_vline(xintercept = fthr, lty = 2) +
    labs(x = "log2(fold change)", y = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 24, color = "black"),
          legend.position = "none",
          panel.grid = element_blank())
}
