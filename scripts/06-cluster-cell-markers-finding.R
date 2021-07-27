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

enriched_top30 = list() ##saving outputs of the for loop in a list. To do some wrangling
# need to loop. Or one can do lapply too. But need to check the syntax.
for (i in 0:(max(conserved_markers$cluster_id))){
  enriched_top30[[paste0("cluster_",i)]] <- enrichr(top30 %>% select(cluster_id, gene) %>% filter( cluster_id== i) %>%  pull(gene), dbs_Allen_Brain_Atlas_10x_scRNA_2021)
  plotEnrich(enriched_top30[[paste0("cluster_",i)]]$`Allen_Brain_Atlas_10x_scRNA_2021`, showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",title = paste0("top30 ","cluster ",i),xlab = "Allen scRNAseq 2021")
  ggsave(paste0("plots/enrichr_markers_top30/","cluster",i,".pdf"), height = 6, width = 7)
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

enriched_top20 = list() ##saving outputs of the for loop in a list. To do some wrangling
# need to loop. Or one can do lapply too. But need to check the syntax.
for (i in 0:(max(conserved_markers$cluster_id))){
  enriched_top20[[paste0("cluster_",i)]] <- enrichr(top20 %>% select(cluster_id, gene) %>% filter( cluster_id== i) %>%  pull(gene), dbs_Allen_Brain_Atlas_10x_scRNA_2021)
  plotEnrich(enriched_top20[[paste0("cluster_",i)]]$`Allen_Brain_Atlas_10x_scRNA_2021`, showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",title = paste0("top20 ","cluster ",i),xlab = "Allen scRNAseq 2021")
  ggsave(paste0("plots/enrichr_markers_top20/","cluster",i,".pdf"), height = 6, width = 7)
  print(paste0("cluster",i," done"))
}

###NOTE: I checked the pdfs, could assign some clusters to some cell types, but this is not exhaustive. Trying another approach:
##Source: https://bioconductor.org/books/release/OSCA/cell-type-annotation.html#assigning-cell-labels-from-gene-sets
# Tried but complicated...

##Trying this to find the cell types 
# 1. Run the EnrichR for top20 and top30 marker genes with for loop and save them in a list. Then clean the lists and keep only the dataframes for the enrichr outputs with cell types for each clusters.
# 2. Do overlap for cluster by cluster and find which cluster has the same cell types between top20 and top30 enriched lists.
##But this will be complex coding, need to use purrr:: from tidyverse . For now, I will do it manually, side by side...

##OR! Do geneOverlap!
###Testing

library(data.table)

Allen_Brain_Atlas_10x_scRNA_2021 = read.table("data/Allen_Brain_Atlas_10x_scRNA_2021.txt", fill = T, header = F, na.strings = "NA")

Allen_Brain_Atlas_10x_scRNA_2021$V0 = paste(Allen_Brain_Atlas_10x_scRNA_2021$V1,Allen_Brain_Atlas_10x_scRNA_2021$V2,Allen_Brain_Atlas_10x_scRNA_2021$V3,Allen_Brain_Atlas_10x_scRNA_2021$V4,Allen_Brain_Atlas_10x_scRNA_2021$V5,Allen_Brain_Atlas_10x_scRNA_2021$V6,sep ="_")

Allen_Brain_Atlas_10x_scRNA_2021 = Allen_Brain_Atlas_10x_scRNA_2021[,c(983,1:982)]
Allen_Brain_Atlas_10x_scRNA_2021 = Allen_Brain_Atlas_10x_scRNA_2021 %>% select(-V1,-V2,-V3,-V4,-V5,-V6)

rm(tmp)
