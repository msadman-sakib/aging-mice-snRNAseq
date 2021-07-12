library(Seurat)
library(tidyverse)
library(metap)
library(future)

#loading data
annotations = readRDS("Rdata/mouse-gene-annotation.rds")
seurat_integrated_12PC = readRDS("Rdata/seurat_integrated_12PC_res_0.3_Hahn.rds")

# change the current plan to access parallelization
plan("multiprocess", workers = 40)

##setting up RAM
options(future.globals.maxSize = 160000 * 1024^2)

plan()

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated_12PC,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:21), get_conserved)


# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (young_avg_log2FC + old_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

saveRDS(conserved_markers, "Rdata/seurat_integrated_12PC_FindConservedMarkers.rds")
saveRDS(top10, "Rdata/seurat_integrated_12PC_top10ConservedMarkers.rds")