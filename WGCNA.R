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

# transpose the data and prepared dataset for WGCNA analysis 
input_mat = t(expression_normalized_data)

allowWGCNAThreads() # allow multi-threading (optional)
# choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to 20, by = 20))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
  )

