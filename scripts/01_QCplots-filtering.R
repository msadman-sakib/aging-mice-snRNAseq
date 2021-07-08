if (!require("pacman")) install.packages("pacman") 
pacman::p_load(dplyr, Seurat, patchwork,ggplot2,stringr,BiocManager,glmGamPoi,sctransform,rio) #p_load installs all those packages if they are not installed, after updating R.

#library(BiocManager)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(glmGamPoi)
library(sctransform)
library(rio)

# 1. Load the filtered dataset from Dennis, metadata generation and saving into Rdata -----
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


# Now following HBC training material
# Add number of genes per UMI for each cell to metadata
agingnuclei$log10GenesPerUMI <- log10(agingnuclei$nFeature_RNA) / log10(agingnuclei$nCount_RNA)


# Compute percent mito ratio
agingnuclei$mitoRatio <- PercentageFeatureSet(object = agingnuclei, pattern = "^mt-")
agingnuclei$mitoRatio <- agingnuclei@meta.data$mitoRatio / 100


# store mitochondrial percentage in object meta data to use in sctransform.
agingnuclei  <- PercentageFeatureSet(agingnuclei, pattern = "^mt-", col.name = "percent.mt.seurat")

##So basically those two mitochondrial stuffs are similar, just the mitoRatio are divided by 100.


# Create metadata dataframe
metadata <- agingnuclei@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
metadata = metadata %>% unite(conditon_replicate, c(ident.conditions, old.ident), remove = FALSE) ###stupido! Typo!

metadata = metadata %>% rename(condition_replicate = conditon_replicate)

# Drop columns
metadata <- metadata %>% select(-old.ident, -percent.mt.seurat) ##removing percent.mt.seurat as it is basically same as mitoRatio.

# Rename columns
metadata <- metadata %>% 
  rename(nUMI = nCount_RNA,
         nGene = nFeature_RNA,
         sample = ident.conditions)




# Add metadata back to Seurat object
agingnuclei@meta.data <- metadata

# Create .RData object to load at any time
save(agingnuclei, file="Rdata/merged_filtered_seurat_agingMice.RData")



# 2. Assessing quality matrix ----------------

##"We will assess various metrics and then decide on which cells are low quality and should be removed from the analysis:
  
#Cell counts
#UMI counts per cell
#Genes detected per cell
#UMIs vs. genes detected
#Mitochondrial counts ratio
#Novelty"

#Why aren't we checking for doublets? Many workflows use maximum thresholds for UMIs or genes, with the idea that a much higher number of reads or genes detected indicate multiple cells. While this rationale seems to be intuitive, it is not accurate. Also, many of the tools used to detect doublets tend to get rid of cells with intermediate or continuous phenotypes, although they may work well on datasets with very discrete cell types. Scrublet is a popular tool for doublet detection, but we haven't adequately benchmarked it yet. Currently, we recommend not including any thresholds at this point in time. When we have identified markers for each of the clusters, we suggest exploring the markers to determine whether the markers apply to more than one cell type.

#In theory, Dennis already did the filtering of the data, so I don't need to do it. But just for visualization that after filtering, what does the data look like, I will generate these plots...


#1.Visualize the number of cell counts per condition
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

ggsave("plots/QC_cellNumber_condition.pdf",height = 7, width = 5)

#2.Visualize the number of cell counts per sample

metadata %>% 
  ggplot(aes(x=condition_replicate, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells per condition")

ggsave("plots/QC_cellNumber_perSample.pdf",height = 5, width = 7)


#3. Visualize the number UMIs/transcripts per cell per condition
# If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("plots/QC_transcriptPerCellperCondition.pdf",height = 5, width = 7)

#4. Visualize the number UMIs/transcripts per cell per sample

metadata %>% 
  ggplot(aes(group=condition_replicate, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("plots/QC_transcriptPerCellpersample.pdf",height = 5, width = 7)

#5. Genes detected per cell
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
ggsave("plots/QC_genes detected per cell.pdf",height = 5, width = 7)


#6. Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(group=condition_replicate, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
ggsave("plots/QC_genes detected per sample.pdf",height = 5, width = 7)



#7. Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("plots/QC_genes detected per cell boxplot.pdf",height = 5, width = 7)


#8. Visualize the distribution of genes detected per sample via boxplot
metadata %>% 
  ggplot(aes(x=condition_replicate, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("plots/QC_genes detected per sample boxplot.pdf",height = 7, width = 12)


### Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot. Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs.

#9. Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
ggsave("plots/QC_FILTERED_UMIs vs. genes detected.pdf",height = 6, width = 12)


#10.# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
##Not saving it, s mitochondrial cells are already removed. This is redundant plot.


#11. Complexity
#The novelty score is computed by taking the ratio of nGenes over nUMI, to see how complex each cell is. Generally, we expect the novelty score to be above 0.80 for good quality cells.

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave("plots/QC_FILTERED_cellComplexity_novelty score above 0.8.pdf",height = 5, width = 7)


# 3. Filtering Cell level (just here for the reference, I am not doing it as cells ar already filtered) ------

#Cell-level filtering
#Now that we have visualized the various metrics, we can decide on the thresholds to apply which will result in the removal of low quality cells. Often the recommendations mentioned earlier are a rough guideline, and the specific experiment needs to inform the exact thresholds chosen. We will use the following thresholds:
#nUMI > 500
#nGene > 250
#log10GenesPerUMI > 0.8
#mitoRatio < 0.2

# Filter out low quality cells using selected thresholds - these will change with experiment
#DO NOT RUN!
#As dennis already filtered it, I am not gonna filter again, as the QC plots looked okay. 
#filtered_seurat <- subset(x = agingnuclei, 
#                          subset= (nUMI >= 500) & 
#                            (nGene >= 250) & 
#                            (log10GenesPerUMI > 0.80) & 
#                            (mitoRatio < 0.20))


# 4. Filtering Gene level (This is important, to keep only genes which are expressed in 10 or more cells. I will do ) ------
# Extract counts
counts <- GetAssayData(object = agingnuclei, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

##Now, keep only genes which are expressed in 10 or more cells. 

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
agingnuclei_filteredCounts = agingnuclei ##maing a copy of the seurat object

agingnuclei_filteredCounts <- CreateSeuratObject(filtered_counts, meta.data = agingnuclei_filteredCounts@meta.data)

# Save filtered subset to new metadata
metadata_clean <- agingnuclei_filteredCounts@meta.data
metadata_clean = metadata_clean %>% rename(condition_replicate = conditon_replicate) %>% select(-nCount_RNA, -nFeature_RNA)

# make a test plot of Complexity
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

##also boxplot for genes
metadata_clean %>% 
  ggplot(aes(x=condition_replicate, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

##It is basically same with those previous QC plots. so, not generating them again.

# Create .RData object to load at any time
save(agingnuclei_filteredCounts, file="Rdata/seurat_filtered_GeneCounts.RData")

