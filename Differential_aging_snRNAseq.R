library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)

# Load the PBMC dataset
agingnuclei.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
agingnuclei <- CreateSeuratObject(counts = agingnuclei.data, project = "agingnuclei80K", names.field = 2, names.delim = "-",min.cells = 3, min.features = 200) # names.field 2 means, the 2nd element after barcode is the sample name, here 1,2,3 is young, 4,5,6,7 is old.

# Stash cell identity classes
agingnuclei[["old.ident"]] <- Idents(object = agingnuclei)
#check old names
levels(agingnuclei)
# rename to new identity
agingnuclei <- RenameIdents(agingnuclei, '1' = 'young', '2' = 'young', '3' = 'young', '4' = "old", '5' =  "old", "6" =  "old", "7" =  "old")
#check new names
levels(agingnuclei)
# Stash new cell identity classes
agingnuclei[["ident.conditions"]] <- Idents(object = agingnuclei) ###To set young and old conditions in the data frame, inside that seurat object.


agingnuclei.list <- SplitObject(agingnuclei, split.by = "ident.conditions")

agingnuclei.list  <- lapply(X = agingnuclei.list , FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


