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

#############Removing underscore from conditon_replicate in both sce and ei objects###Otherwise the later Regex stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+"))) gives only young/old...not the replicate numbers. also, later joining tables of ei doesn't work
sce$conditon_replicate = sub("_", "",  sce$conditon_replicate)

ei$conditon_replicate = sub("_","",ei$conditon_replicate)



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
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), `[`, 1) ###

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)

# Differential gene expression with DESeq2 -----------------------------------------------------------------------------

# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    conditon_replicate = de_samples)

gg_df <- left_join(gg_df, ei[, c("conditon_replicate", "sample")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, conditon_replicate, sample) 

metadata        

# Subsetting dataset to cluster(s) of interest, as an example, for the first cluster in the list, Astrocytes -----------------------------------
# Generate vector of cluster IDs
clusters <- levels(as.factor(metadata$cluster_id))
clusters
clusters[1]

# Subset the metadata to only the Astrocytes
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$conditon_replicate
head(cluster_metadata)

# Subset the counts to only the Astrocytes
counts <- pb[[clusters[1]]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))

# Create DESeq2 object ----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ sample)

## Principal component analysis
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, intgroup = "conditon_replicate" ) + ggtitle("Astrocytes")
ggsave("DE_analysis_scrnaseq/figures/deseq2-astrocytes-pca.pdf", height = 5, width = 6)

## Hierarchical clustering
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("sample"), drop=F])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

cluster_metadata$sample = as.factor(cluster_metadata$sample)

# Results ---------------------------------------------------------------------------------------------------
# Output results of Wald test for contrast for old vs young 
levels(cluster_metadata$sample)[2] ## Young
levels(cluster_metadata$sample)[1] ## Old

contrast <- c("sample", levels(cluster_metadata$sample)[1], levels(cluster_metadata$sample)[2])

# resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.1) ###Changing to padj <0.1

res <- lfcShrink(dds, 
                 coef =  "sample_young_vs_old",
                 type = "apeglm")

# Table of results for all genes --------------------------------------------------------------------------------

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
# No significant genes were found, atleast in Astrocytes... It might be due to the old samples that converge with the youngs...



##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##Now, continuing for all the clusters....
dir.create("DESeq2")
dir.create("DESeq2/pairwise")

# Function to run DESeq2 and get results for all clusters
## x is index of cluster in clusters vector on which to run function
## A is the sample group to compare
## B is the sample group to compare against (base level)

get_dds_resultsAvsB <- function(x, A, B){
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
  rownames(cluster_metadata) <- cluster_metadata$conditon_replicate
  counts <- pb[[clusters[x]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  #all(rownames(cluster_metadata) == colnames(cluster_counts))        
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ sample)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  
  # Plot PCA
  
  DESeq2::plotPCA(rld, intgroup = "sample")
  ggsave(paste0("results/", clusters[x], "_specific_PCAplot.png"))
  
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Plot heatmap
  png(paste0("results/", clusters[x], "_specific_heatmap.png"))
  pheatmap(rld_cor, annotation = cluster_metadata[, c("sample"), drop=F])
  dev.off()
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  # Plot dispersion estimates
  png(paste0("results/", clusters[x], "_dispersion_plot.png"))
  plotDispEsts(dds)
  dev.off()
  
  # Output results of Wald test for contrast for A vs B
  #contrast_A_B <- c("sample", levels(cluster_metadata$sample)[A], levels(cluster_metadata$sample)[B])
  
  # resultsNames(dds)
  res <- results(dds, 
                 #contrast = contrast_A_B,
                 alpha = 0.1) #chanign to 0.1
  
  res <- lfcShrink(dds, 
                   coef =  "sample_young_vs_old", ###I had to hardcode here to make the code run. Be careful next time. 
                   type = "apeglm")
  # Set thresholds
  padj_cutoff <- 0.1 #Set up cutoff.
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$sample)[A], "_vs_", levels(cluster_metadata$sample)[B], "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$sample)[A], "_vs_", levels(cluster_metadata$sample)[B], "_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  ## ggplot of top genes
  normalized_counts <- counts(dds, 
                              normalized = TRUE)
  
  ## Order results by padj values
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n=20)
  
  
  top20_sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% top20_sig_genes)
  
  gathered_top20_sig <- top20_sig_norm %>%
    gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
  
  gathered_top20_sig <- inner_join(ei[, c("conditon_replicate", "sample" )], gathered_top20_sig, by = c("conditon_replicate" = "samplename"))
  
  ## plot using ggplot2
  ggplot(gathered_top20_sig) +
    geom_point(aes(x = gene, 
                   y = normalized_counts, 
                   color = sample), 
               position=position_jitter(w=0.1,h=0)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("log10 Normalized Counts") +
    ggtitle("Top 20 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$sample)[A], "_vs_", levels(cluster_metadata$sample)[B], "_top20_DE_genes.png"))
  
  if(nrow(sig_res)>=2) { ####add this condition for heatmap and volcano plot
  
          #Heatmap
          # Extract normalized counts for only the significant genes
          sig_norm <- data.frame(normalized_counts) %>%
            rownames_to_column(var = "gene") %>%
            dplyr::filter(gene %in% sig_res$gene)
          
          # Set a color palette
          heat_colors <- brewer.pal(6, "YlOrRd")
          
          # Run pheatmap using the metadata data frame for the annotation
          png(paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$sample)[A], "_vs_", levels(cluster_metadata$sample)[B], "all_genes_heatmap.png"),width = 1200, height = 1000, res = 150)
          pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
                   color = heat_colors, 
                   cluster_rows = T, 
                   show_rownames = F,
                   annotation = cluster_metadata[, c("sample", "cluster_id")], 
                   border_color = NA, 
                   fontsize = 10, 
                   scale = "row", 
                   fontsize_row = 10, 
                   height = 20)        
          dev.off()
          #ggsave(paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$sample)[A], "_vs_", levels(cluster_metadata$sample)[B], "all_genes_heatmap.png"))
          
          #Volcano plot 
          
          ## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
          res_table_thres <- res_tbl %>% 
            mutate(threshold = padj < 0.1 & abs(log2FoldChange) >= 0.1)
          
          ## Volcano plot
          ggplot(res_table_thres) +
            geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
            ggtitle(paste0("Volcano plot of ",clusters[x]," relative to control")) +
            xlab("log2 fold change") + 
            ylab("-log10 adjusted p-value") +
            scale_y_continuous(limits = c(0,50)) +
            theme(legend.position = "none",
                  plot.title = element_text(size = rel(1.5), hjust = 0.5),
                  axis.title = element_text(size = rel(1.25)))                    
          ggsave(paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$sample)[A], "_vs_", levels(cluster_metadata$sample)[B], "all_genes_volcanoplot.png"))
          print(paste0("Analysis for ",clusters[x]," done"))
          print("================================")
          } else print("Heatmap or Volcano plot not made as there is less than 2 significant gene")
}

# Run the script on all clusters comparing Old relative to young
# I tried with 0.05. But now trying with 0.1
map(1:length(clusters), get_dds_resultsAvsB, A = 1, B = 2)















