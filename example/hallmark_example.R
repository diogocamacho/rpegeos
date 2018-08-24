library(tidyverse)
library(tidytext)
library(Matrix)
library(limma)
# library(rpegeos)
library(org.Hs.eg.db)
library(GO.db)


# data ----
load("../data/gse18198.RData")
E <- GSE18198$expression
G <- GSE18198$genes
S <- GSE18198$samples

# differential expression ----
hbp_all_ids <- grep("HPB-ALL", S[, 35])
dmso_ids <- grep("DMSO", S$source_name_ch1)
samh_ids <- grep("SAHM1", S$source_name_ch1)

# limma
expdes <- matrix(0, nrow = ncol(E), ncol = 4) # experimental design

expdes[intersect(hbp_all_ids, dmso_ids), 1] <- 1
expdes[intersect(hbp_all_ids, samh_ids), 2] <- 1
expdes[setdiff(dmso_ids, intersect(hbp_all_ids, dmso_ids)), 3] <- 1
expdes[setdiff(samh_ids, intersect(hbp_all_ids, samh_ids)), 4] <- 1
colnames(expdes) <- c("hbp_dmso", "hbp_sahm", "kopt_dmso", "kopt_sahm")

contmat <- makeContrasts(hbp_sahm - hbp_dmso, levels = expdes) # contrast matrix

limma_fit <- lmFit(object = E, design = expdes)
limma_fit <- contrasts.fit(fit = limma_fit, contrasts = contmat)
limma_fit <- eBayes(fit = limma_fit)

res <- topTable(fit = limma_fit,
                coef = 1,
                number = nrow(E),
                sort.by = "none",
                adjust.method = "fdr") # dge results

# volcano plot
res %>%
  ggplot() +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val)), color = "gray") +
  geom_point(data = subset(res, adj.P.Val < 0.05 & abs(logFC) > 1), aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
  labs(x = "log2(Fold change)", y = "-log10(q-value)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20, face = "bold"))

# hallmark pathways ----
load("data/hallmark_pathways.RData")
hallmark_names <- hallmark$names
hallmark_matrix <- hallmark$matrix

# tf-idf ----
tfidf_dtm <- rpegeos_tfidf(data_matrix = hallmark_matrix)
cpm <- crossprod_matrix(tfidf_matrix = tfidf_dtm)


# match genes to pathway
qgs <- G$ENTREZID[which(res$logFC > 1 & res$adj.P.Val < 0.05)] # <--- query gene set
genes_per_pathway <- pathway_genes(gene_set = qgs, pathway_tfidf = tfidf_dtm)

nix1 <- which(genes_per_pathway$number_genes == 0)
if (length(nix1) != 0) {
  ntf <- tfidf_dtm[-nix1, ]
}

nix2 <- which(colSums(ntf) == 0)
if (length(nix2) != 0) {
  ntf <- ntf[, -nix2]
}
ncp <- crossprod_matrix(tfidf_matrix = ntf)

# enrichment ----
# query set
query_vector <- integer(length = ncol(ntf))
query_vector[which(colnames(ntf) %in% qgs)] <- 1

# compute query similarities
query_similarities <- cosine_similarity(tfidf_matrix = ntf,
                                        query_vector = query_vector,
                                        tfidf_crossprod_mat = ncp)

# generate random set
rset <- random_sets(number_sets = 10000,
                    universe_size = ncol(ntf),
                    gene_set_size = sum(query_vector))

# random similarities
asis <- proc.time()
rnd_res <- random_similarity(random_set = rset,
                             tfidf_matrix = ntf,
                             tfidf_crossprod_mat = ncp)
asis <- proc.time() - asis

# dc <- proc.time()
# xxx3 <- do.call(what = cosine_similarity, args = c(query_vector = rset[1:10], list(tfidf_matrix = tfidf_dtm, tfidf_crossprod_mat = cpm)))
# dc <- proc.time() - dc


# compute significance
prandom <- sig_sim(query_similarities, rnd_res)

# global score
rpegeos_score <- (query_similarities / max(query_similarities))  - log10(prandom) + 1
rpegeos_score[prandom == 0] <- (query_similarities[prandom == 0] / max(query_similarities)) - log10(1 / ncol(rnd_res)) + 1


# results
res_enr <- data_frame(pathway_name = hallmark_names[-nix1],
                      cosine_similarity = as.vector(query_similarities),
                      probability_random = prandom,
                      score = as.vector(rpegeos_score),
                      number_genes = genes_per_pathway$number_genes[-nix1],
                      genes_matched = genes_per_pathway$gene_names[-nix1])

# res_enr <- res_enr %>% tibble::add_column(., pathway_name = hallmark_names, .before = 1)

# add pathway names and pathway source
# if (ncol(pathway_names) == 1) {
#   res <- res %>% tibble::add_column(., pathway_name = pathway_names, .before = 1)
# } else if (ncol(pathway_names) == 2) {
#   res <- res %>%
#     tibble::add_column(., pathway_source = pathway_names[,2], .before = 1) %>%
#     tibble::add_column(., pathway_name = pathway_names[,1], .before = 2)
# }
#
# res <- enrichment_results(query_similarities = query_similarities,
#                           random_probability = prandom,
#                           number_genes = genes_per_pathway$number_genes,
#                           gene_names = genes_per_pathway$gene_names,
#                           pathway_names = hallmark_names)


res_enr <- enrich_geneset(gene_set = qgs,
                                  pathway_names = hallmark_names,
                                  tfidf_matrix = tfidf_dtm,
                                  tfidf_crossproduct = cpm,
                                  num_random = 10000)

res_enr %>%
  dplyr::filter(., probability_random < 0.001) %>%
  dplyr::arrange(., desc(cosine_similarity)) %>%
  dplyr::slice(1:15)

plot_enrichment(results_df = res_enr, random_thr = 0.001, max_pathways = 15)
