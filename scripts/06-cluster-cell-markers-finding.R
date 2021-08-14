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

# In Allen paper, they used Top 50. So, trying for Top50
# Now trying with Top 50 genes.
# Extract top 50 markers per cluster
top50 <- conserved_markers %>% 
  mutate(avg_fc = (young_avg_log2FC + old_avg_log2FC) /2) %>%
  group_by(cluster_id) %>%
  top_n(n = 50,
        wt = avg_fc)

enriched_top50 = list() ##saving outputs of the for loop in a list. To do some wrangling
# need to loop. Or one can do lapply too. But need to check the syntax.
for (i in 0:(max(conserved_markers$cluster_id))){
  enriched_top50[[paste0("cluster_",i)]] <- enrichr(top50 %>% select(cluster_id, gene) %>% filter( cluster_id== i) %>%  pull(gene), dbs_Allen_Brain_Atlas_10x_scRNA_2021)
#  plotEnrich(enriched_top50[[paste0("cluster_",i)]]$`Allen_Brain_Atlas_10x_scRNA_2021`, showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",title = paste0("top50 ","cluster ",i),xlab = "Allen scRNAseq 2021")
#  ggsave(paste0("plots/enrichr_markers_top50/","cluster",i,".pdf"), height = 6, width = 7)
  print(paste0("cluster",i," done"))
}

###Trying to remove the overlapped celltypes that are "Human" and "down"
###Testing for manipulating one element.
enriched_top50$cluster_0$Allen_Brain_Atlas_10x_scRNA_2021 %>% filter(!str_detect(Term, "Human")) %>% filter(!str_detect(Term, "down")) %>% view()

enriched_top50_human_down_removed = list()
for (i in 1:(max(conserved_markers$cluster_id)+1)){
  enriched_top50_human_down_removed[[paste0("cluster_",i)]] <- enriched_top50[[paste0("cluster_",i)]]$Allen_Brain_Atlas_10x_scRNA_2021 %>% filter(!str_detect(Term, "Human")) %>% filter(!str_detect(Term, "down"))
    plotEnrich(as.data.frame(enriched_top50_human_down_removed[[i]]), showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",title = paste0("top50 human/down removed ","cluster ",i),xlab = "Allen scRNAseq 2021")
    ggsave(paste0("plots/enriched_top50_human_down_removed/","cluster",i,".pdf"), height = 6, width = 7)
  print(paste0("cluster",i," done"))
}

## It works! Now I have all plots with only the up genes from Allex 10x data! 










###NOTE: I checked the pdfs, could assign some clusters to some cell types, but this is not exhaustive. Trying another approach:
##Source: https://bioconductor.org/books/release/OSCA/cell-type-annotation.html#assigning-cell-labels-from-gene-sets
# Tried but complicated...

##Trying this to find the cell types 
# 1. Run the EnrichR for top20 and top30 marker genes with for loop and save them in a list. Then clean the lists and keep only the dataframes for the enrichr outputs with cell types for each clusters.
# 2. Do overlap for cluster by cluster and find which cluster has the same cell types between top20 and top30 enriched lists.
##But this will be complex coding, need to use purrr:: from tidyverse . For now, I will do it manually, side by side...

##OR! Do geneOverlap!

#####################Testing############

Allen_Brain_Atlas_10x_scRNA_2021 = import("data/Allen_Brain_Atlas_10x_scRNA_2021.csv", fill = T, header = F, na.strings = "NA") ##I had to take the .txt file, put it to google sheets and save as csv, then this import function works.
rownames(Allen_Brain_Atlas_10x_scRNA_2021) = Allen_Brain_Atlas_10x_scRNA_2021$V1
Allen_Brain_Atlas_10x_scRNA_2021 = Allen_Brain_Atlas_10x_scRNA_2021 %>% select(-V1,-V2)

#dataframe to lists
tmp = as.data.frame(t(Allen_Brain_Atlas_10x_scRNA_2021))
allen.genelist =  list()      

for(i in 1:ncol(tmp)) {             # Using for-loop to add columns to list
  allen.genelist[[i]] <- tmp[ , i]
}
names(allen.genelist) <- colnames(tmp)
allen.genelist = lapply(allen.genelist, function(z){ z[!is.na(z) & z != ""]})

#seperate mouse and human marker genes
allen.genelist.mouse = names(allen.genelist) %>% 
                            str_detect('Mouse') %>%
                            keep(allen.genelist, .)

allen.genelist.human = names(allen.genelist) %>% 
  str_detect('Human') %>%
  keep(allen.genelist, .)
#####But it might take long to solve it. Trying other solutions. 

######## Therefore, I will just check the pdfs for Top30 and assign then labels for the obvious ones manually..... 
# So, I did it, I took the top significant term for each cluster from Enrichr(top30). Made an excel file. Here it is.
# Major NOTE!!!!!:::
# So, I took the top significant from Enrichr(using top50 genes per cluster), then used the supplemental file from the paper (https://www.sciencedirect.com/science/article/pii/S0092867421005018, it is in data/Allen_ref_clusters_10.1016-j.cell.2021.04.021.xlsx) and grabbed the cluster details. For now, I focused on the types of cells like gabaergic, glutametargic etc...But there were something weird. I saved the note in a text file in data folder. read it if needed!!

#Then made the excel file that I am loading down below...

clusterID.allen = import("data/manual-top50-celltype-cluster.xlsx")

#load seurat integrated dataset. 
seurat_integrated = readRDS("Rdata/seurat_integrated_PC40_res0.3.rds")

# Rename all identities
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0"="0 Glutamatergic DG 1",
                                  "1"="1 Glutamatergic CA1-do",
                                  "2"="2 Oligo",
                                  "3"="3 GABAergic Sncg 1",
                                  "4"="4 GABAergic Ntng1 HPF 1",
                                  "5"="5 Glutamatergic SUB",
                                  "6"="6 GABAergic Pax6",
                                  "7"="7 Glutamatergic CA3-do",
                                  "8"="8 Oligo-OPC 1",
                                  "9"="9 GABAergic Sncg 2",
                                  "10"="10 Glutamatergic DG 2",
                                  "11"="11 Glutamatergic NP SUB",
                                  "12"="12 Astro",
                                  "13"="13 Glutamatergic DG 3",
                                  "14"="14 Glutamatergic L2 IT ENTl",
                                  "15"="15 GABAergic Ntng1 HPF 2",
                                  "16"="16 Glutamatergic Car3",
                                  "17"="17 Oligo-Glutamatergic",
                                  "18"="18 Micro 1",
                                  "19"="19 GABAergic Lamp5",
                                  "20"="20 GABAergic Sst",
                                  "21"="21 Glutamatergic Mossy",
                                  "22"="22 Glutamatergic L2 IT APr",
                                  "23"="23 Glutamatergic CR",
                                  "24"="24 Gabaergic Glutamatergic",
                                  "25"="25 GABAergic Sst Lamp5 Lhx6",
                                  "26"="26 Peri Endo",
                                  "27"="27 Micro 2",
                                  "28"="28 Micro 3",
                                  "29"="29 Micro 4",
                                  "30"="30 Oligo 3",
                                  "31"="31 Oligo-OPC 2")

# Plot the UMAP
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 2,
        repel = TRUE)
ggsave("plots/final-plots/UMAP_clusters_relabelled.pdf", height = 6, width = 12)

saveRDS(seurat_integrated, "Rdata/seurat_integrated_cell_labelled.rds")


