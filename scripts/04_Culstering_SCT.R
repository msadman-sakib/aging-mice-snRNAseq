#Goals:
  
#To generate cell type-specific clusters and use known cell type marker genes to determine the identities of the clusters.
#To determine whether clusters represent true cell types or cluster due to biological or technical variation, such as clusters of cells in the S phase of the cell cycle, clusters of specific batches, or cells with high mitochondrial content.

#Challenges:
  
#Identifying poor quality clusters that may be due to uninteresting biological or technical variation
#Identifying the cell types of each cluster
#Maintaining patience as this can be a highly iterative process between clustering and marker identification (sometimes even going back to the QC filtering)

#Recommendations:
  
#Have a good idea of your expectations for the cell types to be present prior to performing the clustering. Know whether you expect cell types of low complexity or higher mitochondrial content AND whether the cells are differentiating
#If you have more than one condition, it's often helpful to perform integration to align the cells
#Regress out number of UMIs (by default with sctransform), mitochondrial content, and cell cycle, if needed and appropriate for experiment, so not to drive clustering
#Identify any junk clusters for removal or re-visit QC filtering. Possible junk clusters could include those with high mitochondrial content and low UMIs/genes. If comprised of a lot of cells, then may be helpful to go back to QC to filter out, then re-integrate/cluster.
#If not detecting all cell types as separate clusters, try changing the resolution or the number of PCs used for clustering


# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(ggthemes)

seurat_integrated <- readRDS("Rdata/integrated_seurat_2.5hrs.rds")

# First, need to identify significant PCs(But if using SCTransform, it is not important. I still keep this chunk for reference) -------
## Two methods ------
### 1. Explore heatmap of PCs --------
DimHeatmap(seurat_integrated, 
           dims = 1:12, 
           cells = 500, #The cells argument specifies the number of cells with the most negative or postive PCA scores to use for the plotting. 
           balanced = TRUE)
#The idea is that we are looking for a PC where the heatmap starts to look more "fuzzy", i.e. where the distinctions between the groups of genes is not so distinct.
#This method can be slow and hard to visualize individual genes if we would like to explore a large number of PCs. In the same vein and to explore a large number of PCs, we could print out the top 10 (or more) positive and negative genes by PCA scores driving the PCs.

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

### 2. Elbow plot --------
# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 50)

##########################################
#While the above 2 methods were used a lot more with older methods from Seurat for normalization and identification of variable genes, they are no longer as important as they used to be. This is because the SCTransform method is more accurate than older methods.
##########################################
##########################################
#In theory, with SCTransform, the more PCs we choose the more variation is accounted for when performing the clustering, but it takes a lot longer to perform the clustering. Therefore for this analysis, we will use the first 40 PCs to generate the clusters.
##########################################


# Clustering cells ------
#We will use the FindClusters() function to perform the graph-based clustering. The resolution is an important argument that sets the "granularity" of the downstream clustering and will need to be optimized for every individual experiment. For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.

#The FindClusters() function allows us to enter a series of resolutions and will calculate the "granularity" of the clustering. This is very helpful for testing which resolution works for moving forward without having to run the function for each resolution.

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

tested.resolutions <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.4,1.6,1.8,2.0,4,6,8,10,12,14,16,18,20)

# Determine the clusters for various resolutions. Need produce plots for all of them.
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = tested.resolutions)

# Explore resolutions, this are added as columns in the meta.data.
seurat_integrated@meta.data %>% 
  View()

##So, I have 23 resolutions. How to plot them all together in a grid? I will just make a bar plot to see how many clusters each of those resolutions find. 
#Also, in the meta data, theres a column called seurat cluster. It finds 31 clusters. What is the meaning of it? 
resolutions_clusters =  select(seurat_integrated@meta.data, contains("integrated_snn_res")) 
resolutions_clusters[] = lapply(resolutions_clusters, as.numeric)
resolutions_clusters.summary=data.frame(apply(resolutions_clusters,2,max))
colnames(resolutions_clusters.summary) = "total_clusters"
resolutions_clusters.summary$resolution = rownames(resolutions_clusters.summary)
resolutions_clusters.summary$resolution = gsub("integrated_snn_res.","",as.character(resolutions_clusters.summary$resolution))
resolutions_clusters.summary = resolutions_clusters.summary %>% arrange(as.numeric(resolution))
resolutions_clusters.summary$resolution <- factor(resolutions_clusters.summary$resolution, levels = resolutions_clusters.summary$resolution)

#Line plot
ggplot(data = resolutions_clusters.summary, aes(x = resolution, y=total_clusters, group = 1)) + 
  geom_line(size =1,linetype = "dashed" ) + 
  geom_point(size=2) + scale_y_continuous(name="No. of clusters", breaks = seq(0, 220, by = 10)) + 
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 25) ) #+
  #theme_linedraw()
ggsave("plots/resolutions_no_of_clusters.pdf", height = 8, width = 12)

##Based on this plot, 0.3 and 0.6 are the least resolutions with stable clusters. Let's make cluster UMAP plot for these two only!!

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.3"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3)
ggsave("plots/resolutions_0.3_32_clusters.pdf", height = 5, width = 7)

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3)

ggsave("plots/resolutions_0.6_38_clusters.pdf", height = 5, width = 7)

##But feeling, 0.2 might be better
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.2"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3)

ggsave("plots/resolutions_0.2_27_clusters.pdf", height = 5, width = 7)

