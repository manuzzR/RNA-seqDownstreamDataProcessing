############################################################
##################### Install packages #####################
############################################################

repos <- c("https://cran.rstudio.com/", 
           "https://bioconductor.org/packages/3.16/bioc")
packages = c("SummarizedExperiment", "Biobase", "DESeq2", "annotate",
             "GenomicRanges", "ggdendro", "ggplot2", "grid", "RColorBrewer",
             "vegan", "heatmap3", "lattice", "edgeR", "pheatmap", "reshape2",
             "tidyverse", "svglite")
for (i in packages) {
  if (!i %in% rownames(installed.packages())) {
    install.packages(i, repos = repos)
  }
}
for (i in packages) library(i, character.only = T)

############################################################
#####################  Create dataset  #####################
############################################################

samples <- c('Cp_light3A', 'Cp_light3B', 'Cp_light3C',
            'Cp_Dark3A', 'Cp_Dark3C',
            'Cp_light4A', 'Cp_light4B', 'Cp_light4C',
            'Cp_Dark4A', 'Cp_Dark4B', 'Cp_DarK4C')

conditions <- data.frame(lighting = c(rep ('Light', 3), rep ('Dark', 2),
                                     rep ('Light', 3), rep ('Dark', 3)),
                        concentration = factor(c(rep ('2mg', 5), rep ('4mg', 6))))
conditions$concentration <- relevel(conditions$concentration, ref = "4mg")

# Change the path if needed
countMatrix <- read.delim("DESeq_data_For_Cp_data.txt", header = T, row.names = 1)

rownames(conditions) <- samples
names(countMatrix) <- samples

# The following line filters the count matrix: for each gene at least 2 samples need to have read count > 3
final_counts <- countMatrix[rowSums(countMatrix > 3) >= 2, ]

# In the design, the last condition is never shown, but you can get it in the results.
d.deseq <- DESeqDataSetFromMatrix(countData = final_counts, 
                                  colData = conditions,
                                  design = ~ 0 + lighting + concentration + lighting:concentration)

############################################################
##################  Exploratory analysis  ##################
############################################################

# The argument parallel = parallel::detectCores() will allow it to run faster
dds <- DESeq(d.deseq)  

vsd <- varianceStabilizingTransformation(dds)
theme_set(theme_bw()) # removes the gray background of the following graphs
plotPCA(vsd, intgroup = colnames(conditions))

############################################################
#########################  Results  ########################
############################################################

##all result ##
all_DEG <- results(dds) %>% as.data.frame()
head(all_DEG)

all_DEG %>% select(log2FoldChange,padj) %>%
  filter(row.names(all_DEG) =="TRINITY_DN2246_c1_g1")


# Comparing light vs dark
resultLighting <- results(dds, contrast = c("lighting", "Light", "Dark"))





# This filters the result (the %>% is a pipe, uses the result from the left as input to the right)
# Here we take only genes with absolute log2 fold change > 2 and adjusted p-Value < 0.01
LightingDEGs <- resultLighting %>% as.data.frame() %>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

dim(LightingDEGs)


nameOfUpregulatedGeneInLightVsDark <- row.names(LightingDEGs %>%
                                                  filter(log2FoldChange > 2))

nameofDownregulatingGeneInLightVsDark <- row.names(LightingDEGs %>% 
                                                     filter(log2FoldChange < -2))

factor.labeling.upregulated.lightVsDark <- DataFrame(condition = "upregulated_LightVsDark",
                             gene =nameOfUpregulatedGeneInLightVsDark )

factor.labeling.downregulated.lightVsDark <- DataFrame(condition = "downregulated_LightVsDark",
                                                     gene = nameofDownregulatingGeneInLightVsDark )
                                                                                                          
write.table(x = , file = "factorLabeling.txt", 
            append = F, sep = "\t", quote = F, col.names = F, row.names = F)

write.table( factor.labeling.upregulated.lightVsDark, file = "factor.labeling.upregulated.lightVsDark.txt", 
            append = F, sep = "\t", quote = F, col.names = F, row.names = F)

write.table( factor.labeling.downregulated.lightVsDark, 
             file = "factor.labeling.downregulated.lightVsDark.txt", 
             append = F, sep = "\t", quote = F, col.names = F, row.names = F)



### name list of gene in upregulated and down regulatd gene on light 2mg Vs 4 mg ######

resultConcentration <- results(dds, contrast = c("concentration", "2mg", "4mg"))
ConcentrationDEGs <- resultConcentration %>% as.data.frame() %>% 
  filter(abs(log2FoldChange) > 1 & padj < 0.01)

# For this one, the contrast is the list, so 2mg/L vs 4mg/L in light 
### name list of gene in upregulated and down regulatd gene on light 2mg Vs 4 mg ######
resultLight2mg <- results(dds, list("concentration2mg", "lightingLight.concentration2mg"))
Light2mgDEGs <- resultLight2mg %>% as.data.frame() %>% 
  filter(abs(log2FoldChange) > 2 & padj < 0.01)

### name list of upregulated gene in Ligh2mgDEGs ####

nameOfUpregulatedGeneInLight2mgDEG <- row.names(Light2mgDEGs %>%
                                                  filter(log2FoldChange > 2))
nameOfDownregulatedGeneInLight2mgDEG <- row.names(Light2mgDEGs %>% 
                                                  filter(log2FoldChange < 2))

factor.labeling.upregulated.2mgVs4mg <- DataFrame(condition = "upregulated_2mgVs4mgLight",
                             gene = nameOfUpregulatedGeneInLight2mgDEG)
factor.labeling.Downregulatd.2mgVs4mg <- DataFrame(condition = "downregulatd_2mgvs4mgLight",
                                                  gene = nameOfDownregulatedGeneInLight2mgDEG)

write.table(x =factor.labeling.upregulated.2mgVs4mg , 
            file = "1UpGeneInLight2mgDEGfactorLabeling.txt", 
            append = F, sep = "\t", quote = F, col.names = F, row.names = F)
write.table(x =factor.labeling.Downregulatd.2mgVs4mg, 
            file = "2DownGeneInLight2mgDEGfactorLabeling.txt", 
            append = F, sep = "\t", quote = F, col.names = F, row.names = F)

### genarating diagram from Goseq result, #####
#after Goseq i copy GOseq enriched result on my working directory ###

library(readxl)
library(dplyr)
library(ggplot2)

plot_GOseq <- function(data, number){ as.data.frame(read_tsv (data)) %>%
    na.omit() %>% 
    dplyr::select(category, over_represented_pvalue, 
                  under_represented_pvalue,numDEInCat,numInCat,
                  term ,ontology, over_represented_FDR) %>%
    dplyr::filter(term != 'biological_process') %>%
    dplyr::filter(term != 'cellular_component') %>%
    dplyr::filter(term != 'molecular_function') %>%
    top_n(number, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc,
               y=term,
               size=numDEInCat)) + geom_point() + expand_limits(x=0) +
    labs(x="Hits (%)", y="GO term", colour="p value", size="Count", tittle = data) + 
    theme_bw()+ 
    ggtitle(data)
}




plot_GOseq(data = "downregulated_LightVsDark.GOseq.enriched", number = 20)
plot_GOseq(data = "upregulated_LightVsDark.GOseq.enriched", number = 20)




