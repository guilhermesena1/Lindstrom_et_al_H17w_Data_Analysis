# Lindstrom et al Human Fetal Kidney Data Analysis

This repository reproduces the figures for the paper "Progressive recruitment of mesenchymal progenitors highlights a time-dependent process of cell fate acquisition in mouse and humannephrogenesis"

The scripts that generated the figures are available in the file H17w_Data_Analysis.Rmd and can be run on Rstudio.

Please note that results may vary slightly because no random seed was set for the figures that went into the paper and some methods have random steps to it (Seurat Clustering, Jackstraw and GMM). However, our reruns of the methods have shown that the results are consistent despite potential changes in the number of nephron and interstitial progenitor clusters.  

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
Datasets are available in the 10x sparse matrix format in GEO
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112570
(download the file GSE112570_RAW.tar	in the bottom of the page)

Sparse matrices can either be read using the `readMM` function in the Matrix package or the `Read10X` function in Seurat.

If you wish to reanalyze the data from the raw reads (eg: Map it to a different reference or run quality control metrics), you can download the .bam files with tagged cell barcodes and UMIs directly from SRA:

https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6921770 (Kidney 1)
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6921771 (Kidney 2)

These can be converted back to fastq to run on upstream processing tools (see the [10x guidelines](https://support.10xgenomics.com/docs/bamtofastq) for instructions)
