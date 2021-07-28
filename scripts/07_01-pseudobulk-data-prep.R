# July 2021

### psuedobulk data prep
# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)


# Bring in Seurat object
#seurat = readRDS("Rdata/seurat_integrated_cell_labelled.rds")
seurat = seurat_integrated

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts

metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident) ###This is to get the cluster identities, I can set it to whatever I want!!

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "orig.ident")] ##sample_id is orig.ident here.

##Let's create directories for pseudobulk DE analysis
# dir.create("DE_analysis_scrnaseq")
# dir.create("DE_analysis_scrnaseq/data")
# dir.create("DE_analysis_scrnaseq/results")
# dir.create("DE_analysis_scrnaseq/figures")

saveRDS(sce, "DE_analysis_scrnaseq/data/scRNA-seq_input_data_for_DE.rds")
