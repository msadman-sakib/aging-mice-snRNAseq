## Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
library(Seurat)
library(cowplot)
library(tidyverse)
library(metap)
library(rio)
library(enrichR)
#load data #previously was using seurat_integrated_12PC_res_0.3_Hahn.rds, but now I decided to use seurat_integrated_PC40_res0.3.rds after checking the line plot and UMAP embedding. 

seurat_integrated = readRDS("Rdata/seurat_integrated_PC40_res0.3.rds")

n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample" )) %>%
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

n_cells_proportion = n_cells %>%  pivot_longer(!sample, names_to = "clusters", values_to = "count")
n_cells_proportion$clusters = factor(n_cells_proportion$clusters, levels = c(0:(max(as.numeric(seurat_integrated@meta.data$integrated_snn_res.0.3))-1))) ##changing ggplot x axis order
old = n_cells_proportion %>% filter(sample == "old") %>% mutate(old_norm = count/4)
young = n_cells_proportion %>% filter(sample == "young") %>% mutate(old_norm = count/3)
n_cells_proportion = full_join(old, young)
n_cells_proportion =n_cells_proportion %>% rename(norm.count = old_norm)
ggplot(n_cells_proportion, aes(x = clusters, y = norm.count, fill = sample)) +
  geom_col(position = "fill") + geom_hline(yintercept=0.5, linetype="dashed", color = "black", size = 1)
ggsave("plots/cluster_proportions-Normed-young-old.pdf", width = 10, height = 7) ##But I have also generated unnormed proportions. Check that as well.


# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
ggsave("plots/cluster_proportions-young-old-UMAP.pdf", width = 14, height = 7)

# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()


# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
ggsave("plots/cluster_QC_metrics.pdf", height = 15, width = 15) ##VERY IMPORTANT PLOT!!!

#(disregard the cluster numbers here)The metrics seem to be relatively even across the clusters, with the exception of the nUMIs and nGene exhibiting higher values in clusters 3, 9, 14, and 15, and, perhaps, cluster 17. We will keep an eye on these clusters to see whether the cell types may explain the increase.If we see differences corresponding to any of these metrics at this point in time, then we will often note them and then decide after identifying the cell type identities whether to take any further action.

###############Exploration of the PCs driving the different clusters
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)


# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
ggsave("plots/cluster_QC_16PCs.pdf", height = 15, width = 15) ##VERY IMPORTANT PLOT!!!
#(disregard the cluster numbers here)We can see how the clusters are represented by the different PCs. For instance, the genes driving PC_2 exhibit higher expression in clusters 6, 11, and 17 (maybe a bit higher in 15, too). We could look back at our genes driving this PC to get an idea of what the cell types might be:

# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)
#####Real answer!! PC_1 positive genes here are Plp1, Pde4b, Mbp, St18, Prr5l . Mbp is oligodendrocyte gene marker. So, cluster 2 should be oligodendrocyte!!! yay! First cluster identified!!

##########To truly determine the identity of the clusters and whether the resolution is appropriate, it is helpful to explore a handful of known gene markers for the cell types expected.


#Exploring known cell type markers
#NOTE: The SCTransform normalization was performed only on the 3000 most variable genes, so many of our genes of interest may not be present in this data.

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

##from CoolMPS paper (Mature Oligodendrocytes = Plp1, Mog, Pde4b, St18, Slc24a2, Pcdh9; choroid plexus cells = Ttr, Htr2c; Astrocytes = Apoe, Slc1a2, Slc1a3, Prex2, Cst3, Gabrb1, Gpc5; Neurons = Lsamp, Kcnip4, Tenm2, Nrg3, Celf2, Grin2a).

#Depending on our markers of interest, they could be positive or negative markers for a particular cell type. The combined expression of our chosen handful of markers should give us an idea on whether a cluster corresponds to that particular cell type.

#For the markers used here, we are looking for positive markers and consistency of expression of the markers across the clusters. For example, if there are two markers for a cell type and only one of them is expressed in a cluster - then we cannot reliably assign that cluster to the cell type.

##From naim vai: 16 July 2021 Marker genes.....Try it..
#Neuron: Syt1, Rbfox3, Scl17a7, Gad2
#Astrocytes: Gja1, Ndrg2
#Oligodendrocyte: Mal, Cldn11,Mbp
#Microglia: Dock8, Ikzf1, Hexb, Inpp5d


#Mature Oligo markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("Plp1", "Mog", "Pde4b", "St18", "Slc24a2", "Pcdh9"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
ggsave("plots/cluster_QC-40PC-mature-olig.pdf", width = 12, height = 12)


#Astro markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("Apoe", "Slc1a2", "Slc1a3", "Prex2", "Cst3", "Gabrb1", "Gpc5"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
ggsave("plots/cluster_QC-40PC-Astro.pdf", width = 12, height = 12)

#Microglia markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("Tmem119", "Itgam", "Ptprc","Aif1","Cx3cr1","Cd68","Cd40"), #Itgam Cd11b was not found
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
ggsave("plots/cluster_QC-40PC-microglia.pdf", width = 12, height = 12)

####NOTE: If any cluster appears to contain two separate cell types, it's helpful to increase our clustering resolution to properly subset the clusters. Alternatively, if we still can't separate out the clusters using increased resolution, then it's possible that we had used too few principal components such that we are just not separating out these cell types of interest. To inform our choice of PCs, we could look at our PC gene expression overlapping the UMAP plots and determine whether our cell populations are separating by the PCs included.

###########I ran this chuck on HPC using slurm. Script name: 05-01-Findconservedmarker-HPC
# 
# 
# #It is tough like that. I could only identify oligodendrocytes. Therefore, need to Identify conserved markers in all conditions. As I have two conditions, need to run FindConservedMarkers
# DefaultAssay(seurat_integrated_12PC) <- "RNA"
# 
# ##Lots of reading to do again in 09_merged_SC_marker_identification.md
# 
# #loading data
# annotations = readRDS("Rdata/mouse-gene-annotation.rds")
# seurat_integrated_12PC = readRDS("Rdata/seurat_integrated_12PC_res_0.3_Hahn.rds")
# 
# # Create function to get conserved markers for any given cluster
# get_conserved <- function(cluster){
#   FindConservedMarkers(seurat_integrated_12PC,
#                        ident.1 = cluster,
#                        grouping.var = "sample",
#                        only.pos = TRUE) %>%
#     rownames_to_column(var = "gene") %>%
#     left_join(y = unique(annotations[, c("gene_name", "description")]),
#               by = c("gene" = "gene_name")) %>%
#     cbind(cluster_id = cluster, .)
# }
# 
# # Iterate function across desired clusters
# system.time(conserved_markers <- map_dfr(c(0:21), get_conserved)) 
# 
# 
# # Extract top 10 markers per cluster
# top10 <- conserved_markers %>% 
#   mutate(avg_fc = (young_avg_log2FC + old_avg_log2FC) /2) %>% 
#   group_by(cluster_id) %>% 
#   top_n(n = 10, 
#         wt = avg_fc)
# 
# # Visualize top 10 markers per cluster
# View(top10)
# 
# 
# saveRDS(conserved_markers, "Rdata/seurat_integrated_12PC_ConservedMarkers.rds")
# saveRDS(top10, "Rdata/seurat_integrated_12PC_top10ConservedMarkers.rds")

#Finding markers for all clusters
#it would take quite a while to run. Also, it is possible that when you run this function on all clusters, in some cases you will have clusters that do not have enough cells for a particular group - and your function will fail. For these clusters you will need to use FindAllMarkers().
###########I ran this chuck on HPC. Script name: sakib-test-conserved.R

#Now, loading the outputs from HPC outputs

conserved_markers = readRDS("Rdata/seurat_integrated_12PC_FindConservedMarkers.rds")
top10 = readRDS("Rdata/seurat_integrated_12PC_top10ConservedMarkers.rds")

export(top10, "results/top10ConservedMarkers.xlsx")
export(conserved_markers, "results/allConservedMarkers.xlsx")


#Extract top 20 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (young_avg_log2FC + old_avg_log2FC) /2) %>%
  group_by(cluster_id) %>%
  top_n(n = 20,
        wt = avg_fc)

export(top20, "results/top20ConservedMarkers.xlsx")

#Extract top 30 markers per cluster
top30 <- conserved_markers %>% 
  mutate(avg_fc = (young_avg_log2FC + old_avg_log2FC) /2) %>%
  group_by(cluster_id) %>%
  top_n(n = 30,
        wt = avg_fc)

export(top30, "results/top30ConservedMarkers.xlsx")

###Use Toppfun/Toppcell atlas to see!!

##Or use EnrichR!

#I am using EnrichR and using Allen scRNAseq atlas data in a different script, script 06.



##Now, if this doesn't work(say, some cluster not giving any genes) then one needs to use FindMarker. But I have already atleast 30 significant genes in each clusters.