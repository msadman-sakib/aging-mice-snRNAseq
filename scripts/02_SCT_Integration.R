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

#The following like of code took around 15 mins
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
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

# Integration using CCA ----------
# Load the split seurat object into the environment if needed
#split_seurat <- readRDS("data/split_seurat.rds")

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)


# Now, we are going to perform CCA, find the best buddies or anchors and filter incorrect anchors. This might take half an hour.

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

## The goal of integration is to ensure that the cell types of one condition/dataset align with the same celltypes of the other conditions/datasets (e.g. control macrophages align with stimulated macrophages).

















# Splitting young and old
agingnuclei.list <- SplitObject(agingnuclei, split.by = "ident.conditions")


##Normalization by sctransform, using "glmGamPoi" method to make calculations faster. ----
agingnuclei.list  <- lapply(X = agingnuclei.list , FUN = function(x) {
  x <- SCTransform(x, method = "glmGamPoi", verbose = T, vars.to.regress = "percent.mt")
})

##pre-requisites before integration/anchoring
features <- SelectIntegrationFeatures(object.list = agingnuclei.list, nfeatures = 3000)
agingnuclei.list <- PrepSCTIntegration(object.list = agingnuclei.list, anchor.features = features)
#save.image("TempTilsctransformNorm.Rdata")

## so far this worked -----

##for rpca, need to run pca individuall first
agingnuclei.list <- lapply(X = agingnuclei.list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = F)
})

##---sofar this worked---#### took 1-2 min

# finding anchors cell by cell between young and old ------
agingnuclei.anchors <- FindIntegrationAnchors(object.list = agingnuclei.list, normalization.method = "SCT", anchor.features = features,reduction = "rpca", dims = 1:50) ##This takes forever with the default pca. use rpca for less time. rpca worked nicely! seems like it is working faster than normal "pca" as reduction argument.total time only 2 mins! Also, if the integration doesnt overlap well in UMAP, one can change " k.anchor = 20" in FindIntegrationAnchors function to make them overlapping.
##---sofar this worked---#### took 2 min

agingnuclei.combined.sct <- IntegrateData(anchorset = agingnuclei.anchors, normalization.method = "SCT", verbose = T, dims = 1:50)##adding dims argument in both here and last line seemed to make it work

###comment...so here, due to making the analysis fast, I did the rpca method according to here: https://satijalab.org/seurat/articles/integration_large_datasets.html
#I only followed the rpca method, did not use the reference method though. In theory according to them, CCA can be superior if, "....the datasets are highly divergent (for example, cross-modality mapping or cross-species mapping), where only a small subset of features can be used to facilitate integration, and you may observe superior results using CCA..."..

##---sofar this worked---#### took 5 minutes.

# this plot is important to determine number of PCs to include in dims flag for clustering...
ElbowPlot(agingnuclei.combined.sct)

##clustering
agingnuclei.combined.sct <- RunPCA(agingnuclei.combined.sct, verbose = T)
agingnuclei.combined.sct <- RunUMAP(agingnuclei.combined.sct, reduction = "pca", dims = 1:50)

###need to do these two step for cell clustering numbers, resolution controls number of clusters.
agingnuclei.combined.sct <- FindNeighbors(agingnuclei.combined.sct, reduction = "pca", dims = 1:50)
agingnuclei.combined.sct <- FindClusters(agingnuclei.combined.sct, resolution = 0.5)

###plotting
p1 <- DimPlot(agingnuclei.combined.sct, reduction = "umap", group.by = "ident.conditions")
p2 <- DimPlot(agingnuclei.combined.sct, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
ggsave("plots/young_old_cell_clusters.pdf", height = 5, width = 10)

##---sofar this worked!! got the UMAP plot for old vs young and cell clusters numbers... next, find markers and diff expression
save.image("TempTillUMAPplots.Rdata")


# next, identifying conserved cell type markers ------

# This only got markers among default 3000 variable genes -------
# For performing differential expression after integration, we switch back to the original
# data
#DefaultAssay(agingnuclei.combined.sct) <- "RNA" ##can be also "integrated", it contains the transformed values..But RNA was used in this tutorial: https://satijalab.org/seurat/archive/v3.1/immune_alignment.html.
#But it was not changed here in this updated tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#finding-differentially-expressed-features-cluster-biomarkers-  ##following this 2nd one..

# find markers for every cluster compared to all remaining cells, report only the positive ones -------

######Took 1.5 hours for this next line!!!
agingnuclei.markers <- FindAllMarkers(agingnuclei.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) ##It took 1.5 hours!!!

agingnuclei.markers.top2 = agingnuclei.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) ##seeing top 2 genes in each markers.

save(x=agingnuclei.markers.top2, file="agingnuclei.markers.table.Rdata")

agingnuclei.markers.top = agingnuclei.markers %>% group_by(cluster)
save(x=agingnuclei.markers.top, file="agingnuclei.markers.tableAll.Rdata")


##Now inspecting top genes in each cluster
top.markers= agingnuclei.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC)

###test
FeaturePlot(agingnuclei.combined.sct, features = top.markers %>% pull(gene) %>% .[1:6] , min.cutoff = "q9")

FeaturePlot(agingnuclei.combined.sct, features = c("Mbp"), min.cutoff = "q1")

##But this is not correct, as it is based on only 3000 genes. I need all genes. Here, for example, I searched for Rbfox3, and it was not there! Therefore, I need to change the DefaultAssay() to "RNA" or "SCT"....

# Finding marker correctly with all genes  ------

DefaultAssay(agingnuclei.combined.sct) <- "RNA"
###$!!! Probably need to normalize first, then run the FindAllMarker here!!!!
agingnuclei.markers.RNA <- FindAllMarkers(agingnuclei.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #Now it is showing longer minutes for each cluster...around 2.5hours!
save(x=agingnuclei.markers.RNA, file="agingnuclei.markers.RNA.took3hours.Rdata")


agingnuclei.markers.RNA.top2 = agingnuclei.markers.RNA %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) ##seeing top 2 genes in each markers.
export(agingnuclei.markers.RNA.top2, "agingnuclei.markers.top2.xlsx")

agingnuclei.markers.RNA.topall = agingnuclei.markers.RNA %>% group_by(cluster)
export(agingnuclei.markers.RNA.topall, "agingnuclei.markers.tableAll.xlsx")

#testing some gene expressions plots, now it has information for 20K genes. So my genes are there.
FeaturePlot(agingnuclei.combined.sct, features = c("Mog", "Rbfox3", "Gfap", "Icam1"), min.cutoff = "q1")

plots <- VlnPlot(agingnuclei.combined.sct, features = c("Mog", "Rbfox3", "Gfap", "Icam1"), split.by = "ident.conditions", group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

##these are okay for few genes. But need a heatmap for more genes
##since in defaultassay = RNA, the scale.data is empty, need to do it.
agingnuclei.combined.sct = NormalizeData(agingnuclei.combined.sct) ##this is happening in the RNA assay.
all.genes <- rownames(agingnuclei.combined.sct)
agingnuclei.combined.sct <- ScaleData(agingnuclei.combined.sct, features = all.genes) 

agingnuclei.markers.RNA %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top3 ##here, I made Top1, top5 manually and generated the plots. But top 3 heatmap looks good. 

pdf("plots/top3_byclusterMarker.pdf", width = 27, height = 10)
DoHeatmap(agingnuclei.combined.sct, features = top3$gene) + NoLegend()
dev.off()

#save.image("TemptillClusterPlots.Rdata")


#### but now using sctransform data to find the marker genes. For this, making a new assay for SCT for all genes. Lets see. If not works, then can remove by setting the assay NULL. 
##But I just read that, doing this will fod sure crash as it needs huge memors...one possibility is to run on cluster.

# STILL UNDER DEVELOPMENT, about using SCTransform to scale.data for all genes, then use it for UMAP plot...Maybe I will not need it! -------


agingnuclei.combined.sct.marker.genes = SCTransform(agingnuclei.combined.sct, method = "glmGamPoi",assay = "RNA",new.assay.name = "SCTallGenes",return.only.var.genes = F, verbose = T, vars.to.regress = "percent.mt")

#set to new assay
DefaultAssay(agingnuclei.combined.sct.marker.genes) <- "SCTallGenes"

#find marker genes
agingnuclei.markers.SCT <- FindAllMarkers(agingnuclei.combined.sct.marker.genes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot ="scale.data") #Now it is showing longer minutes for each cluster...around 2.5hours!

agingnuclei.markers.SCT %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top3

DoHeatmap(agingnuclei.combined.sct, features = top3$gene) + NoLegend()

FeaturePlot(agingnuclei.combined.sct, features = c("Mbp"), min.cutoff = "q1")

# Cell type annotation, then subsetting neurons and glia for diff analysis -----
##But now it is quite difficult to annotate each clusters. 
# need to use AddModuleScore() function to add the cell type marker genes in the agingnuclei.combined.sct object. 

#example:
#agingnuclei.combined.sct <- AddModuleScore(agingnuclei.combined.sct,
#                             features = list(vector with genes),
#                            name="give a name that will be used to plot the list, like neurons, OPC, microglia etc...")

ElbowPlot(agingnuclei.combined.sct, ndims = 50) #to see the PCs for selecting for clustering...

##Date: 7 July 2021
#I am pausing this script for today. I need to restart the analysis comprehensively to do marker identification and cell type labelling. Need major update of the code. That's why, will start a new script under the same repo. 