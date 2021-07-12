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
library(clustree)

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



#So, from previous script, 03_01_checking_cell_clusters_to_remove_sparse_clusters.R, I will use 40 pca for UMAP projection, as there is no difference between 40 and 50, and even reducing to 20 or 30, the small clusters stays, meaning, they are real!!



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


##But feeling, 0.1 might be better
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.1"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3)

ggsave("plots/resolutions_0.1_18_clusters.pdf", height = 5, width = 7)



##But this method is not exhaustive. I found clustree package to deal with this elegantly.

library(clustree)
##check the prefix of the different resolutions

seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = tested.resolutions)


head(seurat_integrated[[]])

#making clustree
clustree(seurat_integrated, prefix = "integrated_snn_res.")
ggsave("plots/resolutions_clustree_23res.pdf", height = 10, width = 20)
##But the plot looks super hazy. Need to simplify. I need to use less resolutions.

#setting res from 0.1-1.0, as they have less than 50 clusters. 
seurat_integrated.clustree <- seurat_integrated ##making a copy

##removing those extra res to exclude from clustree

seurat_integrated.clustree@meta.data$integrated_snn_res.1.4 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.1.6 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.1.8 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.2 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.4 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.6 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.8 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.10 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.12 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.14 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.16 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.18 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.20 = NULL

clustree(seurat_integrated.clustree, prefix = "integrated_snn_res.")
ggsave("plots/resolutions_clustree_10res.pdf", height = 10, width = 20)

##now removing other res to make clustree plot from 0.1-0.6.
seurat_integrated.clustree@meta.data$integrated_snn_res.0.8 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.0.9 = NULL
seurat_integrated.clustree@meta.data$integrated_snn_res.1 = NULL
clustree(seurat_integrated.clustree, prefix = "integrated_snn_res.")
ggsave("plots/resolutions_clustree_10res.pdf", height = 10, width = 20)


################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
######################################################
#But still, cannot decide. Now I checked Oliver Hahn paper, they used FindNeighbors fimd 1:12. In my elbow plot, at 15th PC, the Standard Dev is the lowest, then the change is not that much. This is similar to 12 PC used in oliver hahn paper. So, now generating using 15 PCs.

##In the mean time, I tried for 15PC- see elbow plot...But the result is still not convincing. I think I will follow Oliver hahn values, 12 PC, 0.4 Res...

# Determine the K-nearest neighbor graph
seurat_integrated_12PC <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:12)

tested.resolutions_12PC <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

# Determine the clusters for various resolutions. Need produce plots for all of them.
seurat_integrated_12PC <- FindClusters(object = seurat_integrated_12PC,
                                  resolution = tested.resolutions_12PC)


seurat_integrated_12PC@meta.data$integrated_snn_res.1.4 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.1.6 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.1.8 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.2 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.4 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.6 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.8 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.10 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.12 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.14 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.16 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.18 = NULL
seurat_integrated_12PC@meta.data$integrated_snn_res.20 = NULL

##So, I have 10 resolutions. How to plot them all together in a grid? I will just make a bar plot to see how many clusters each of those resolutions find. 
#Also, in the meta data, theres a column called seurat cluster. It finds 31 clusters. What is the meaning of it? 
resolutions_clusters_12PC =  select(seurat_integrated_12PC@meta.data, contains("integrated_snn_res")) 
resolutions_clusters_12PC[] = lapply(resolutions_clusters_12PC, as.numeric)
resolutions_clusters_12PC.summary=data.frame(apply(resolutions_clusters_12PC,2,max))
colnames(resolutions_clusters_12PC.summary) = "total_clusters"
resolutions_clusters_12PC.summary$resolution = rownames(resolutions_clusters_12PC.summary)
resolutions_clusters_12PC.summary$resolution = gsub("integrated_snn_res.","",as.character(resolutions_clusters_12PC.summary$resolution))
resolutions_clusters_12PC.summary = resolutions_clusters_12PC.summary %>% arrange(as.numeric(resolution))
resolutions_clusters_12PC.summary$resolution <- factor(resolutions_clusters_12PC.summary$resolution, levels = resolutions_clusters_12PC.summary$resolution)

#Line plot
ggplot(data = resolutions_clusters_12PC.summary, aes(x = resolution, y=total_clusters, group = 1)) + 
  geom_line(size =1,linetype = "dashed" ) + 
  geom_point(size=2) + scale_y_continuous(name="No. of clusters", breaks = seq(0, 220, by = 10)) + 
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 25) ) #+
#theme_linedraw()
ggsave("plots/resolutions_12PC_no_of_clusters.pdf", height = 8, width = 12)



# Assign identity of clusters
Idents(object = seurat_integrated_12PC) <- "integrated_snn_res.0.4"
# Plot the UMAP
DimPlot(seurat_integrated_12PC,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

ggsave("plots/12PC_resolutions_0.4_like Oliver_hahn.pdf", height = 5, width = 7)

##Trying 0.3
# Assign identity of clusters
Idents(object = seurat_integrated_12PC) <- "integrated_snn_res.0.3"
# Plot the UMAP
DimPlot(seurat_integrated_12PC,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

ggsave("plots/12PC_resolutions_0.3_like Oliver_hahn.pdf", height = 5, width = 7)

##Trying 0.2
# Assign identity of clusters
Idents(object = seurat_integrated_12PC) <- "integrated_snn_res.0.2"
# Plot the UMAP
DimPlot(seurat_integrated_12PC,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

ggsave("plots/12PC_resolutions_0.2_like Oliver_hahn.pdf", height = 5, width = 7)


##With this, I got 22 clusters, same as Oliver Hahn...So, I will keep this values.
###12 PC for cell embeding, 0.3 for cluster identification...

#making clustree
clustree(seurat_integrated_12PC, prefix = "integrated_snn_res.")
ggsave("plots/12PC_resolutions_clustree_10res.pdf", height = 10, width = 20)
##but honestly, it is really not easy to use to fix which res to use...saving as a reference...


# Finally, assign identity of clusters to 0.3, for 12 PC cell embedding, NOT 40 PC embedding!!
Idents(object = seurat_integrated_12PC) <- "integrated_snn_res.0.3"

# Plot the UMAP
DimPlot(seurat_integrated_12PC,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
ggsave("plots/12PC_resolutions_0.3_like Oliver_hahn_FINAL.pdf", height = 8, width = 12)

saveRDS(seurat_integrated_12PC, file = "Rdata/seurat_integrated_12PC_res_0.3_Hahn.rds")

##IMPORTANT Note: 10July, 01:04 AM
#I tried today to selct the PCs for embedding and resolutions for cluster number selection. I finally fixated using the oliver Hahn parameters from CoolMPS paper. But now, I realised, I should be following the tutorial for this. THey have a QC step for this. Damn! Still, I learned a lot today! Selecting PC and resolution is very crucial for downstream analysis, therefore, it is good that I tried all these combinations. In the next script, I will do the clustering QC. If that works, than Alhamdulillah! If not, then I will come back here, change PCs and resolution probably. 
#But also realized, this embedding(FindNeighbors(), ndims=, PC selection) and resolution(FindClusters(), resolution= ) has to be just precise, cause there is no right or wrong. Once I check via the marker genes, it can be that neurons are dispersed across cluster 1,2,3, so those 3 clusters can be named as neurons, for example. 