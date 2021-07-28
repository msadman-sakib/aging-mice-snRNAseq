# Single-cell RNA-seq - normalization

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

#load data
load(file = "Rdata/seurat_filtered_GeneCounts.RData")

#Normalization and regressing out unwanted variation --------
  #Goals:
  #To accurately normalize and scale the gene expression values to account for differences in sequencing depth and overdispersed count values.
  #To identify the most variant genes likely to be indicative of the different cell types present.

  #Recommendations:
  #Have a good idea of your expectations for the cell types to be present prior to performing the clustering. Know whether you expect cell types of low complexity or higher mitochondrial content AND whether the cells are differentiating
  #Regress out number of UMIs (default using sctransform), mitochondrial content, and cell cycle, if needed and appropriate for experiment, so not to drive clustering downstream

# 1.Seeing cell cycle effects(No effect) -------

#The most common biological data correction is to remove the effects of the cell cycle on the transcriptome. 

#First we roughly normalize the data, we do SCTransform later...This normalization is solely for the purpose of exploring the sources of variation in our data.

# Normalize the counts
seurat_phase <- NormalizeData(agingnuclei_filteredCounts)

#load mouse cell cycle marker genes
load("Rdata/cell_cycle_markers_Mmu.Rdata")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)                                


#After scoring the cells for cell cycle, we would like to determine whether cell cycle is a major source of variation in our dataset using PCA. To perform PCA, we need to first choose the most variable features, then scale the data. Since highly expressed genes exhibit the highest amount of variation and we don't want our 'highly variable genes' only to reflect high expression, we need to scale the data to scale variation with expression level. The Seurat ScaleData() function will scale the data by:

#adjusting the expression of each gene to give a mean expression across cells to be 0
#scaling expression of each gene to give a variance across cells to be 1

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE) ##these are the defaults for FindVariableFeatures()

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)


# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")


# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
ggsave("plots/QC_cell-cylce-gene-expression-no-effect.pdf",height = 6, width = 10)

#We do not see large differences due to cell cycle phase. Based on this plot, we would not regress out the variation due to cell cycle.

# 2.Seeing mitochondrial expression effects(some effect) ------
#    Check quartile values
summary(seurat_phase@meta.data$mitoRatio)
#But all the quartiles are Zeros, as I already removed the cells with high Mitochondria values. So, I can make a cut at Mean.

# Turn mitoRatio into categorical factor vector based on Mean value
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 5.024e-05, Inf), 
                                     labels=c("Nonexistent","very low(>0.005024%)"))


# Plot the PCA colored by mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr",
        )
ggsave("plots/QC_Mitochondria-gene-expression-some-effect.pdf",height = 6, width = 10)

#There is a clear difference. So we need to regress out the mitochondrial expression


# 3.Regressing out unwanted variables using SCTransform(in this case, Mitochondria) -----

# Sctransform automatically accounts for cellular sequencing depth by regressing out sequencing depth (nUMIs). However, if there are other sources of uninteresting variation identified in the data during the exploration steps we can also include these. We observed little to no effect due to cell cycle phase and so we chose not to regress this out of our data. We observed some effect of mitochondrial expression and so we choose to regress this out from the data.

#Since we have two samples in our dataset (from two conditions), we want to keep them as separate objects and transform them as that is what is required for integration. We will first split the cells in seurat_phase object into "Young" and "Old":

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

split_seurat <- split_seurat[c("young", "old")]

#Before we run this for loop, we know that the output can generate large R objects/variables in terms of memory. If we have a large dataset, then we might need to adjust the limit for allowable object sizes within R (Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:
options(future.globals.maxSize = 8000 * 1024^2) ##setting 8GB.

# For loop for doing SCTransfrom for all samples.

#NOTE: By default, after normalizing, adjusting the variance, and regressing out uninteresting sources of variation, SCTransform will rank the genes by residual variance and output the 3000 most variant genes. If the dataset has larger cell numbers, then it may be beneficial to adjust this parameter higher using the variable.features.n argument.

!!!!!!!!!!!!!!!###########NOTE!!!


#The following like of code took around 15 mins
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio")) ####NOTE: Should also include ‘nCount_RNA’, ‘nFeature_RNA’ !!!!!!
}
##Note: This can be done using lapply and function too. See Seurat Vignette. 
##I had warnings, more than 50. All had this same line
#1: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
##But dont worry. It is a non-issue!! From here: https://github.com/ChristophH/sctransform/issues/25
#"These warnings are showing that there are some genes for which it is hard to reliably estimate theta (presumably because of very few non-zero observations). Usually we don't worry about these warnings too much, since we regularize the parameters in a later step, thus averaging out uncertainty of individual gene parameters"


# Check which assays are stored in objects
split_seurat$young@assays

# Save the split seurat object
saveRDS(split_seurat, "Rdata/split_seurat.rds")
