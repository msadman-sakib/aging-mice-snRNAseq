# Single-cell RNA-seq - normalization

# Integration using CCA ----------
## The goal of integration is to ensure that the cell types of one condition/dataset align with the same celltypes of the other conditions/datasets (e.g. control macrophages align with stimulated macrophages).

# Load the split seurat object into the environment if needed
#split_seurat <- readRDS("data/split_seurat.rds")
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)

# change the current plan to access parallelization----ONLY FOR THE CLUSTER

options(future.globals.maxSize= 160000*1024^2)
print("RAM usage set to 160 GB)")

# Load the split seurat object into the environment if needed
split_seurat <- readRDS("split_seurat.rds")
print("split_seurat.rds loaded")

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
print("SelectIntegrationFeatures Ran successfully")

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

print("PrepSCTIntegration Ran successfully")

# Now, we are going to perform CCA, find the best buddies or anchors and filter incorrect anchors. This might take half an hour-1hour.
#########################################################
print("Find best buddies - can take a while to run, Running FindIntegrationAnchors")
plan("multiprocess", workers = 40) ##ONLY FOR RUNNING IN CLUSTER
print("multiprocess set to 40 CPUs)")
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
print("FindIntegrationAnchors Finished!")


# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
print("IntegrateData Finished!")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
print("RunPCA Finished!")

# Plot PCA
PCAPlot(seurat_integrated,
        group.by = "sample")  
ggsave("Integrated-PCA.pdf", height = 7, width = 14)
print("PCA plot saved!")

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
print("RunUMAP  Finished!")

# Plot UMAP                             
DimPlot(seurat_integrated, group.by = "sample",)     
ggsave("Integrated-UMAP.pdf", height = 7, width = 7)
print("UMAP plot saved!")

# Plot UMAP split by sample
DimPlot(seurat_integrated,
        group.by = "sample",
        split.by = "sample")
ggsave("Integrated-UMAP-sidebyside.pdf", height = 7, width = 14) ###here just for a dot, I have to rerun the whole thing shit!!!
print("Integrated side by side UMAP plot saved!")

# Save integrated seurat object
saveRDS(seurat_integrated, "integrated_seurat.rds")
print("integrated_seurat.rds saved")

