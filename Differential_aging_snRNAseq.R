if (!require("pacman")) install.packages("pacman") 
pacman::p_load(dplyr, Seurat, patchwork,ggplot2,stringr,BiocManager ) #p_load installs all those packages if they are not installed, after updating R.

#library(BiocManager)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(glmGamPoi)
library(sctransform)


# Load the dataset
agingnuclei.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
agingnuclei <- CreateSeuratObject(counts = agingnuclei.data, project = "agingnuclei80K", names.field = 2, names.delim = "-",min.cells = 3, min.features = 200)
# names.field 2 means, the 2nd element after barcode is the sample name, here 1,2,3 is young, 4,5,6,7 is old.

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

agingnuclei.list <- SplitObject(agingnuclei, split.by = "ident.conditions")

# store mitochondrial percentage in object meta data to use in sctransform.
agingnuclei.list  <- lapply(X = agingnuclei.list , FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^mt-", col.name = "percent.mt")
})

##Normalization by sctransform, using "glmGamPoi" method to make calculations faster. ----
agingnuclei.list  <- lapply(X = agingnuclei.list , FUN = function(x) {
  x <- SCTransform(x, method = "glmGamPoi", verbose = T, vars.to.regress = "percent.mt")
})

save.image("TempTilsctransformNorm.Rdata")
##---sofar this worked---####

##pre-requisites befor integration/anchoring
features <- SelectIntegrationFeatures(object.list = agingnuclei.list, nfeatures = 3000)
agingnuclei.list <- PrepSCTIntegration(object.list = agingnuclei.list, anchor.features = features)

##finding anchors cell by cell
agingnuclei.anchors <- FindIntegrationAnchors(object.list = agingnuclei.list, normalization.method = "SCT", anchor.features = features)
agingnuclei.combined.sct <- IntegrateData(anchorset = agingnuclei.anchors, normalization.method = "SCT")

##clustering
agingnuclei.combined.sct <- RunPCA(agingnuclei.combined.sct, verbose = T)
agingnuclei.combined.sct <- RunUMAP(agingnuclei.combined.sct, reduction = "pca", dims = 1:30)


###plotting
p1 <- DimPlot(agingnuclei.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(agingnuclei.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE, repel = TRUE)
p1 + p2


