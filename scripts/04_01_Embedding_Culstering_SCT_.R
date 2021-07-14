# Load libraries

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(ggthemes)
library(clustree)
library(gridExtra)


##ram setup
options(future.globals.maxSize = 6000 * 1024^2)

#load data
seurat_integrated <- readRDS("Rdata/integrated_seurat_2.5hrs.rds")

#So, from previous script, 03_01_checking_cell_clusters_to_remove_sparse_clusters.R, I will use 40 pca for UMAP projection, as there is no difference between 40 and 50, and even reducing to 20 or 30, the small clusters stays, meaning, they are real!!
##NOTE!!!: 14 July 2021: used 12 PCs and resolution 0.2/0.3 but the embedding was not perfect. It was though to find clusters and markers. That's why, again, need to try 


# Determine the K-nearest neighbor graph - so, embedding and here in dims, I am putting the Cluster amounts. 

#########First, I am doing for 12 PCs
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:12)

tested.resolutions <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

# Determine the clusters for various resolutions. Need produce plots for all of them.
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = tested.resolutions)

# Explore resolutions, this are added as columns in the meta.data.
seurat_integrated@meta.data %>% 
  View()

##So, I have 10 resolutions. How to plot them all together in a grid? I will just make a bar plot to see how many clusters each of those resolutions find. 
#Also, in the meta data, theres a column called seurat cluster. It finds 31 clusters. What is the meaning of it? 

##data.frame making for resolutions and how it effects the number of clusters.
resolutions_clusters =  select(seurat_integrated@meta.data, contains("integrated_snn_res")) 
resolutions_clusters[] = lapply(resolutions_clusters, as.numeric)
resolutions_clusters.summary=data.frame(apply(resolutions_clusters,2,max))
colnames(resolutions_clusters.summary) = "total_clusters"
resolutions_clusters.summary$resolution = rownames(resolutions_clusters.summary)
resolutions_clusters.summary$resolution = gsub("integrated_snn_res.","",as.character(resolutions_clusters.summary$resolution))
resolutions_clusters.summary = resolutions_clusters.summary %>% arrange(as.numeric(resolution))
resolutions_clusters.summary$resolution <- factor(resolutions_clusters.summary$resolution, levels = resolutions_clusters.summary$resolution)

#sample Line plot
p1 = ggplot(data = resolutions_clusters.summary, aes(x = resolution, y=total_clusters, group = 1)) + 
  geom_line(size =1,linetype = "dashed" ) + 
  geom_point(size=2) + scale_y_continuous(name="No. of clusters", breaks = seq(0, 220, by = 10)) + 
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 25) ) #+
#theme_linedraw()
#ggsave("plots/resolutions_12PCs_clusterNumbers.pdf", height = 8, width = 12)

##But after I checked the elbow plot, I think I should also try 30 PCs! But seems, the resolution = 0.3/0.4 should be optimum. Uff! This is really mad! there is no consensus on how to find an appropriate PCs, and resolution!

######################For loop for saving multiple datasets in a list##############

#####################Then do it for 15,30,40 PCS###Can be included for other PCS and resolutions here too! Make this a function and publich it to github!

#load data
seurat_integrated <- readRDS("Rdata/integrated_seurat_2.5hrs.rds")

PCs= c(12,15,30,40) ##set which PCs for embedding
tested.resolutions <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) ##set which resolutions for clustering
#make empty list
PCs.cluster.numbers = list()
for(i in PCs){
  seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                     dims = 1:i)
  
  # Determine the clusters for various resolutions. Need produce plots for all of them.
  seurat_integrated <- FindClusters(object = seurat_integrated,
                                    resolution = tested.resolutions)
  resolutions_clusters =  select(seurat_integrated@meta.data, contains("integrated_snn_res")) 
  resolutions_clusters[] = lapply(resolutions_clusters, as.numeric)
  resolutions_clusters.summary=data.frame(apply(resolutions_clusters,2,max))
  colnames(resolutions_clusters.summary) = paste0("PC_",i)
  resolutions_clusters.summary$resolution = rownames(resolutions_clusters.summary)
  resolutions_clusters.summary$resolution = gsub("integrated_snn_res.","",as.character(resolutions_clusters.summary$resolution))
  resolutions_clusters.summary = resolutions_clusters.summary %>% arrange(as.numeric(resolution))
  resolutions_clusters.summary$resolution <- factor(resolutions_clusters.summary$resolution, levels = resolutions_clusters.summary$resolution)
  
  PCs.cluster.numbers[[paste0("PC_",i)]] = resolutions_clusters.summary
  print(paste0("PC_",i," dataframe saved!"))
}
#using reduce() from tidyverse::purrr to join multiple dataframes into one single dataframe, then gathering the columns with different PCs and the no. of clusters into two columns:
PCs.cluster.numbers.df = PCs.cluster.numbers %>% reduce(left_join, by = "resolution") %>% gather(-resolution, key = "PC_numbers", value = "Number_of_clusters")

#line plot for all those PCs and their corresponding number of clusters, for the defined resolutions. 
ggplot(PCs.cluster.numbers.df) + geom_line(aes(x = resolution, y = Number_of_clusters,group = PC_numbers, colour = PC_numbers)) + scale_y_continuous(breaks = seq(0, 48, by = 2))
ggsave("plots/resolutions-combined-4diff.PCS.pdf", height = 5, width = 6)

######################For loop for saving multiple plots in a list##############

##Okay now, I am gonna use res = 0.3 for PCs, as it will give around/ less than 30 clusters. 
#PCs to test, I am also including here PC=40, which I did in the very beginning!
#delete, then load data, otherwise Rstudio crashes.
#rm(seurat_integrated)
seurat_integrated <- readRDS("Rdata/integrated_seurat_2.5hrs.rds")
PCs = c(12,15,30,40)
PCs.plots = list()
##FOR LOOP TO SAVE MULTIPLE UMAP PLOTS IN A LIST!!!
for(i in PCs ){

  seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:i)
# Determine the clusters for various resolutions. Need produce plots for all of them.
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = 0.3) ##this is fixed to get less than 30 clusters in each PCs.
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:i) ####I MISSED THIS STEP! without this, the FindNeughbours dims won't show effect!!!
PCs.plots[[paste0("PC_",i)]] = DimPlot(seurat_integrated,
                  reduction = "umap",
                  label = TRUE,
                  label.size = 3)
print(paste0("PC_",i," plot saved!"))
}
#now, plotting all those saved in that list in a single grid
g = grid.arrange(grobs =  PCs.plots, ncol = 2)
ggsave(plot = g, "plots/PC_embedding-test-res0.3/all4PCs_new40PC.pdf", height = 12, width = 15)

##Okay now the embedding works, after including the RunUMAP function.

## all4PCs_new40PC.pdf this plot has 4 UMAPS for 4 PCs(12,15,30,40 using 0.3 resolution for clusters. Seems like PC40 is good. But need to decide now...)




###NOTE for troubleshooting: 

#Upps! The embedding did not change! The cells were the same in all three plots. Lets try manually first, then come back here..
#seurat_integrated <- readRDS("Rdata/integrated_seurat_2.5hrs.rds")
# seurat_integrated <- FindNeighbors(object = seurat_integrated, 
#                                    dims = 1:30)
# seurat_integrated <- FindClusters(object = seurat_integrated,
#                                   resolution = 0.3) 
# seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:12)
# DimPlot(seurat_integrated,
#         reduction = "umap",
#         label = TRUE,
#         label.size = 3)
# #dir.create("plots/PC_embedding-test-res0.3")
# ggsave("plots/PC_embedding-test-res0.3/PC30.pdf", height = 5, width = 6)
# rm(seurat_integrated)
# seurat_integrated <- readRDS("Rdata/integrated_seurat_2.5hrs.rds")
# seurat_integrated <- FindNeighbors(object = seurat_integrated, 
#                                    dims = 1:12)
# seurat_integrated <- FindClusters(object = seurat_integrated,
#                                   resolution = 0.3) 
# DimPlot(seurat_integrated,
#         reduction = "umap",
#         label = TRUE,
#         label.size = 3)
# #dir.create("plots/PC_embedding-test-res0.3")
# ggsave("plots/PC_embedding-test-res0.3/PC12.pdf", height = 5, width = 6)
# 
# ##Okay, so big reveal, embedding did not changed. I think I need to reload the data everytime before doing the FindNeughbors. Lets do that. So doing this above..for PC 12..I also forgot to do the RUNUMAP function for embedding! Okay, I found the problem I think!!
# ##Now, running the forloop again, with that RUNUMAP with correct dims! 