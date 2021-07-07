if (!require("pacman")) install.packages("pacman") 
pacman::p_load(dplyr, Seurat, patchwork,ggplot2,stringr,BiocManager,glmGamPoi,sctransform,rio) #p_load installs all those packages if they are not installed, after updating R.

#library(BiocManager)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(glmGamPoi)
library(sctransform)
library(rio)


# Load the dataset
agingnuclei.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
agingnuclei <- CreateSeuratObject(counts = agingnuclei.data, project = "agingnuclei80K", names.field = 2, names.delim = "-",min.cells = 3, min.features = 200) #names.field 2 means, the 2nd element after barcode is the sample name, here 1,2,3 is young, 4,5,6,7 is old.

# Stash cell identity classes
agingnuclei[["old.ident"]] <- Idents(object = agingnuclei)

#check old names
levels(agingnuclei)

# rename to new identity
agingnuclei <- RenameIdents(agingnuclei, '1' = 'young', '2' = 'young', '3' = 'young', '4' = "old", '5' =  "old", "6" =  "old", "7" =  "old")

#check new names
levels(agingnuclei)

# Stash new cell identity classes
agingnuclei[["ident.conditions"]] <- Idents(object = agingnuclei) ###This is to set young and old conditions in the data frame, inside that seurat object.
# Splitting young and old
agingnuclei.list <- SplitObject(agingnuclei, split.by = "ident.conditions")

# store mitochondrial percentage in object meta data to use in sctransform.
agingnuclei.list  <- lapply(X = agingnuclei.list , FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^mt-", col.name = "percent.mt")
})

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

# For performing differential expression after integration, we switch back to the original
# data
#DefaultAssay(agingnuclei.combined.sct) <- "RNA" ##can be also "integrated", it contains the transformed values..But RNA was used in this tutorial: https://satijalab.org/seurat/archive/v3.1/immune_alignment.html.
#But it was not changed here in this updated tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#finding-differentially-expressed-features-cluster-biomarkers-  ##following this 2nd one..

# find markers for every cluster compared to all remaining cells, report only the positive ones -------

######Took 1.5 hours for this next line!!!
DefaultAssay(agingnuclei.combined.sct) <- "RNA"
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

DefaultAssay(agingnuclei.combined.sct) <- "RNA"
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

save.image("TemptillClusterPlots.Rdata")


#### but now using sctransform data to find the marker genes. seems like the 
DefaultAssay(agingnuclei.combined.sct) <- "SCT"
agingnuclei.markers.SCT <- FindAllMarkers(agingnuclei.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #Now it is showing longer minutes for each cluster...around 2.5hours!


# Cell type annotation, then subsetting neurons and glia for diff analysis -----
##But now it is quite difficult to annotate each clusters. 
# need to use AddModuleScore() function to add the cell type marker genes in the agingnuclei.combined.sct object. 

#example:
#agingnuclei.combined.sct <- AddModuleScore(agingnuclei.combined.sct,
#                             features = list(vector with genes),
#                            name="give a name that will be used to plot the list, like neurons, OPC, microglia etc...")


