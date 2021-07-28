library(Seurat)
library(tidyverse)

seurat_integrated = readRDS("Rdata/integrated_seurat_2.5hrs.rds")



# Plot PCA
PCAPlot(seurat_integrated,
        group.by = "sample")  


ElbowPlot(seurat_integrated, ndims= 50)
ggsave("plots/cell-density-Elbowplot.pdf", height = 5, width = 8)


# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40, #default was 1:40
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,group.by = "sample", cols = c("#00BFC4","#F8766D"), order = "old") #need to set colors as the orientation somehow was switched.
ggsave("plots/cell-density_pca40.pdf", height = 5, width = 7)



# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:50, #default was 1:40
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,group.by = "sample", cols = c("#00BFC4","#F8766D"), order = "old") #need to set colors as the orientation somehow was switched.
ggsave("plots/cell-density_pca50.pdf", height = 5, width = 7)



# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:30, #default was 1:40
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,group.by = "sample", cols = c("#00BFC4","#F8766D"), order = "old") #need to set colors as the orientation somehow was switched.
ggsave("plots/cell-density_pca30.pdf", height = 5, width = 7)


# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:20, #default was 1:40
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,group.by = "sample", cols = c("#00BFC4","#F8766D"), order = "old") #need to set colors as the orientation somehow was switched.
ggsave("plots/cell-density_pca20.pdf", height = 5, width = 7)


# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:10, #default was 1:40
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,group.by = "sample", cols = c("#00BFC4","#F8766D"), order = "old") #need to set colors as the orientation somehow was switched.
ggsave("plots/cell-density_pca10.pdf", height = 5, width = 7)



##Like oliver Hahn, to see the cell embedding
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:12, #default was 1:40
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,group.by = "sample", cols = c("#00BFC4","#F8766D"), order = "old") #need to set colors as the orientation somehow was switched.
ggsave("plots/cell-density_pca12.pdf", height = 5, width = 7)



###Stopping here, I will use 40 pca for UMAP projection, as there is no difference between 40 and 50, and even reducing to 20 or 30, the small clusters stays, meaning, they are real!!

###After 3 hours

###Stopping here, I will use 12 pca for cell embedding, to match Oliver Hahn's number of clusters.
