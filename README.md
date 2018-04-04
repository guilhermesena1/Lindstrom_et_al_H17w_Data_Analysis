# Lindstrom et al Human Fetal Kidney Data Analysis

This repository reproduces the figures for the paper "Progressive recruitment of mesenchymal progenitors highlights a time-dependent process of cell fate acquisition in mouse and humannephrogenesis"

The scripts that generated the figures are available in the file H17w_Data_Analysis.Rmd and can be run on Rstudio.

## Dependencies to install
You need 6 packages to run the Rmd.

```
# CRAN Packages
install.packages("Seurat")
install.packages("dplyr")
install.packages("WGCNA")
install.packages("igraph")
install.packages("mclust")

# Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("monocle")
```

## Count table datasets
Datasets are available both in the 10x sparse matrix format in GEO:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112570

We also made available the matrix in the uncompressed tab-separated format at:
xxxxxxxx
