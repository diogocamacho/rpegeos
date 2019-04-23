# rPEGEOS: Pathway Enrichment for Gene ExpressiOn Signatures with R

This package performs gene set enrichment for a given gene expression signature based on the TF-IDF matrix, computed on the pathway representations. rPEGEOS comes pre-packaged with 16 different gene sets on **human** gene sets. Other species can be accomodated. [Contact me](mailto:diogo.camacho@wyss.harvard.edu) and we'll chat. 

## Installing
In R, do:

```r
library(devtools)
devtools::install_github("diogocamacho/rpegeos")
```

## Running rPEGEOS
The easiest way to run rPEGEOS is through its wrapper `enrich_geneset`. All you need is a set of genes that you want to enrich, as Entrez gene IDs.

### Wrapper
To make it easier to run, a wrapper is also deployed.  To enrich a gene set (provided you have a tf-idf matrix and a cross-product vector, do:

```r
ES <- enrich_geneset(gene_set)
```

where `gene_set` is an Nx3 matrix where the columns are:

  - Column 1: Entrez IDs
  - Column 2: Log2(fold change)
  - Column 3: p-value (or adjusted p-value) for differential expression
  
There is no need to pre-filter the gene set, as `rPEGEOS` will do that for you, based on thresholds defined by the user. Follow on-screen instructions and you are all set!


### Plotting et al.
The protocol above will generate all the necessary elements to generate a final results data frame.  The results can be plotted with `ggplot2` or your favorite plotting tools. A plot of the top enrichment results will be generated as a consequence of running the wrapper.


### In the works
An R Shiny package is being developed to perform gene set enrichments as a web service.


