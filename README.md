# rPEGEOS: Pathway Enrichment for Gene ExpressiOn Signatures with R

This package performs pathway enrichment for a given gene expression signature based on tf-idf.


## Installing
In R, do:

```r
library(devtools)
devtools::install_github("diogocamacho/rpegeos")
```

## Running rPEGEOS
rPEGEOS requires the `dplyr`, `Matrix`, and `tidytext` packages from R.  After getting a pathway set collection (for example, from [MSigDB](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) or [CPDB](http://cpdb.molgen.mpg.de)), generate a binary matrix where each row is a pathway and each column is a unique gene. The columns of the pathway-gene matrix will be used by `rPEGEOS` to compute the cosine similarity. A good practice would be to keep identifiers such as the Entrez IDs as the names for the rows. 

### TF-IDF matrix
The tfidf matrix for rPEGEOS can be generated with: 

```r
M <- rpegeos::tfidf(pathway_data)
```

where `pathway_data` is the one-hot encoded matrix.  

### Cross-product
The next step is to calculate the cross-product of the tf-idf matrix, which will be used to determine the cosine similarity between the query vector and the tf-idf matrix. 

```r
cpm <- crossprod_matrix(M)
```


### Query vector
To generate the query vector, do:

```r
query_vector <- integer(length = ncol(tfidf_matrix))
query_vector[which(colnames(tfidf_matrix) %in% gene_set)] <- 1
```

which will generate the one-hot encoded vector.  

### Cosine similarity
Next, we will calculate the cosine similairity between the query vector and each pathway set:

```r
query_similarities <- cosine_similarity(tfidf_matrix = M,
                                            query_vector = query_vector,
                                            tfidf_crossprod_mat = cpm)
```

### Random probabilities
To determine the probability of the cosine similarity between any given pathway set and the query vector being random, do:

```r
prandom <- random_probability(similarity_results = query_similarities,
                                gs_size = sum(query_vector),
                                num_sets = 10000,
                                target_tfidf = M,
                                tfidf_crossprod_mat = cpm)
```

where `num_sets` is the number of random sets to be generated (defaults to 10000).  


### Plotting et al.
The protocol above will generate all the necessary elements to generate a final results data frame.  The results can be plotted with `ggplot2` or your favorite plotting tools.


### Wrapper
To make it easier to run, a wrapper is also deployed.  To enrich a gene set (provided you have a tf-idf matrix and a cross-product vector, do:

```r
ES <- enrich_geneset(gene_set, pathway_names, tfidf_matrix, crossproduct, num_random)
```

where `gene_set` is a vector of gene IDs (again, same ID source as that used for tf-idf matrix -- Entrez, Ensembl, gene symbol... -- just be consistent), `pathway_names` is the names of the pathways in your pathway sets (ie, row names for tf-idf matrix), and the other 3 arguments are as defined above.

