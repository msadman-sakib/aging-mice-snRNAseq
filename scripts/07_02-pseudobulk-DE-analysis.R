# July 2021

### psuedobulk Differential expression analysis with DESeq2
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

sce = readRDS("DE_analysis_scrnaseq/data/scRNA-seq_input_data_for_DE.rds")

# Explore the raw counts for the dataset ------------------------------------------------------------------------------------------------------------------------

## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

counts(sce)[1:6, 1:6] #We don't want to run head() on this dataset, since it will still show the thousands of columns, so we just looked at the first six rows and columns.


## Explore the cellular metadata for the dataset
dim(colData(sce))

head(colData(sce)) #Here conditon_replicate / orig.ident = sample_id, sample = group_id. 

# Acquiring necessary metrics for aggregation across cells in a sample ------------------------------------------------------------------------------------------------------------------------

# Named vector of cell cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$conditon_replicate)))

# Total number of samples 
ns <- length(sids)
ns

## To perform sample-level differential expression analysis, we need to generate sample-level metadata. To do this, we will reorder samples in the single-cell metadata to match the order of the factor levels of the sample ID, then extract only the sample-level information from the first cell corresponding to that sample.

# Generate sample level metadata  ------------------------------------------------------------------------------------------------------------------------

## Determine the number of cells per sample
table(sce$conditon_replicate)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$conditon_replicate))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$conditon_replicate)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei ##This has lots of columns like MitoRatio etc, which I dont need. But for now, not dropping them. 


# Perform QC if not already performed
dim(sce)

# Calculate quality control (QC) metrics
sce <- addPerCellQC(sce, subsets=list(Mito=grep("mt-", rownames(sce)))) #'calculateQCMetrics' is defunct. use addPerCellQC


colnames(colData(sce)) ###New columns are added, called "sum" "detected" "subsets_Mito_sum" "subsets_Mito_detected" "subsets_Mito_percent" "total".
# detected means detected genes, the number of features for the cell that have counts above the detection limit (default of zero)
# sum total number of counts for the cell (i.e., the library size)
# total means? 

####deactivating plot function for now....
# plotColData(sce, x = "sum", y="detected", colour_by="sample") + xlab("total gene counts/cell") + ylab("total detected genes/cell")
# ggsave("DE_analysis_scrnaseq/figures/QC_genes_detectedVsTotalCounts.pdf", height = 6, width = 7)


# Get cells w/ few/many detected genes
sce$is_outlier <- isOutlier(
  metric = sce$detected, ##instead of total_features_by_counts(it meant genes detected per cell). Here, we use detected instead of that. 
  nmads = 2, type = "both", log = TRUE) #nmads: A more flexible way of choosing thresholds is through the isOutlier function. This defines the threshold at a certain number of median absolute deviations (MADs) away from the median. Values beyond this threshold are considered outliers and can be filtered out, assuming that they correspond to low-quality cells. 

##Original cell numbers: 68452
## Using nmads 2: 61018
## Using nmads 3: 68055
## I would say, I stick to using 2. 

# Remove outlier cells
sce <- sce[, !sce$is_outlier]
dim(sce)

## Remove lowly expressed genes which have less than 5 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 5, ]

dim(sce)

##Initial gene count: 22609
##genes which have less than 10 cells with any counts: 13655 
##genes which have less than 5 cells with any counts: 15188 ##I will keep at least 5 cells with any counts for a given gene. 

# Now, we are ready for aggregation of counts to the sample level. Essentially, we are taking the sum of counts for each sample within each cell type.

# Count aggregation to sample level --------------------------------------------------------------------------------------------------------------------------------------------
# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "conditon_replicate")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

pb[1:6, 1:6]

##Now I have individual cells and given replicates!! WOOHOOOOO!!! Alhamdulillah!!

# To perform DE analysis on a per cell type basis, we need to wrangle our data in a couple ways. We need to do the following steps:
#   
# Split our data by cell type
# Transform the matrix so that the genes are the row names and the samples are the column names

# We will split our data by cell type; however, not always do all samples contain cells of every cell type. To determine which samples are present for each cell type we can run the following:
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 4), `[`, 1) 




