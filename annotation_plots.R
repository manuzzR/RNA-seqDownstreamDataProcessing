library(AnnotationDbi)
library(parallel)
library(tidyverse)
library(GO.db)

# Modify number of GOs to show
topN = 10
#

goannotation = read.delim("C:/analysis_with_karen/go_annotations.txt", header = F)
names(goannotation) = c("Genes", "GO_terms")
allGOterms =
  goannotation$GO_terms %>%
  strsplit(split = ",") %>%
  unlist() %>%
  unname %>%
  unique()

tableGOTerms = mclapply(
  allGOterms, mc.cores = detectCores(), function(x){
  grepl(x, goannotation$GO_terms) %>% sum
})

tableGOs = data.frame(GO = allGOterms[1:10], 
                      occurences = unname(unlist(tableGOTerms))
                      )
tableGOs$Ontology = GOTERM[tableGOs$GO] %>% Ontology
tableGOs$Term = GOTERM[tableGOs$GO] %>% Term


tables_top = list()
for (i in c("BP", "MF", "CC")) {
  tables_top[[i]] = filter(tableGOs, Ontology == i)
  tables_top[[i]] = tables_top[[i]][ order(tables_top[[i]]$occurences, 
                                           decreasing = T), ]
  showGON = ifelse(nrow(tables_top[[i]]) < topN, nrow(tables_top[[i]]), topN)
  tables_top[[i]] = tables_top[[i]][ 1:showGON, ]
  tables_top[[i]]$Term = as.factor(tables_top[[i]]$Term)
}

tableGOplot = rbind(tables_top[[1]], tables_top[[2]], tables_top[[3]]) %>%
  data.frame

write.delim(tableGOplot, x = "tableGOPlot.txt", row.names = F, 
            col.names = T, append = F, quote = F, sep = "\t")