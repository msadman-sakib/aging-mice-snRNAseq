##Need to manually check what genes corresponds to what celltypes...It needs lot of prior knowledge on marker genes..
!!!!!!!!!!()###Last would be using AddModule score and check gene list expression...

#Now, loading the outputs from HPC outputs
library(Seurat)
library(tidyverse)
library(rio)
library(enrichR)
#load dataset
conserved_markers = readRDS("Rdata/seurat_integrated_FindConservedMarkers.rds")

# Testing cluste cell type idenfication manually using enrichr. -------

# Enrichr settings 
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs %>%  view() ##list of databases to choose from.
#setting up the database. This contains the Allen brain single cell Atlas gene types.
dbs_Allen_Brain_Atlas_10x_scRNA_2021 <- c("Allen_Brain_Atlas_10x_scRNA_2021")

# Extract top 30 markers per cluster 
top30 <- conserved_markers %>% 
  mutate(avg_fc = (young_avg_log2FC + old_avg_log2FC) /2) %>%
  group_by(cluster_id) %>%
  top_n(n = 30,
        wt = avg_fc)
# need to loop. Or one can do lapply too. But need to check the syntax.
for (i in 0:(max(conserved_markers$cluster_id))){
  enriched_top30 <- enrichr(top30 %>% select(cluster_id, gene) %>% filter( cluster_id== i) %>%  pull(gene), dbs_Allen_Brain_Atlas_10x_scRNA_2021)
  
  plotEnrich(enriched_top30$`Allen_Brain_Atlas_10x_scRNA_2021`, showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value")
  
  ggsave(paste0("plots/enrichr_markers_top30/","cluster",i,".pdf"))
  print(paste0("cluster",i," done"))
}

#Some clusters are obvious. But not all of them. Lets try with top 20 and top 10 genes. 
#dir.create("plots/enrichr_markers_top20/")

# Now trying with Top 20 genes.
# Extract top 20 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (young_avg_log2FC + old_avg_log2FC) /2) %>%
  group_by(cluster_id) %>%
  top_n(n = 20,
        wt = avg_fc)
# need to loop. Or one can do lapply too. But need to check the syntax.
for (i in 0:(max(conserved_markers$cluster_id))){
  enriched_top20 <- enrichr(top20 %>% select(cluster_id, gene) %>% filter( cluster_id== i) %>%  pull(gene), dbs_Allen_Brain_Atlas_10x_scRNA_2021)
  
  plotEnrich(enriched_top20$`Allen_Brain_Atlas_10x_scRNA_2021`, showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value")
  
  ggsave(paste0("plots/enrichr_markers_top20/","cluster",i,".pdf"))
  print(paste0("cluster",i," done"))
}

###NOTE: I checked the pdfs, could assign some clusters to some cell types, but this is not exhaustive. Trying another approach:
##Source: https://bioconductor.org/books/release/OSCA/cell-type-annotation.html#assigning-cell-labels-from-gene-sets
# Tried but complicated...
