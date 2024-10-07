# uploading nessesary packages

repos <- c("https://cran.rstudio.com/", 
           "https://bioconductor.org/packages/3.16/bioc")
packages = c("SummarizedExperiment", "Biobase", "DESeq2", "annotate",
             "GenomicRanges", "ggdendro", "ggplot2", "grid", "RColorBrewer",
             "vegan", "heatmap3", "lattice", "edgeR", "pheatmap", "reshape2",
             "tidyverse", "svglite", "wgcna")
for (i in packages) {
  if (!i %in% rownames(installed.packages())) {
    install.packages(i, repos = repos)
  }
}
for (i in packages) library(i, character.only = T)
