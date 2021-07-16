library(Seurat)
library(tidyverse)
library(metap)
library(future)

#loading data
annotations = readRDS("Rdata/mouse-gene-annotation.rds")
seurat_integrated = readRDS("Rdata/seurat_integrated_PC40_res0.3.rds")

# change the current plan to access parallelization
plan("multiprocess", workers = 40)

##setting up RAM
options(future.globals.maxSize = 160000 * 1024^2)

plan()

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:(max(as.numeric(seurat_integrated@meta.data$integrated_snn_res.0.3))-1)), get_conserved)
saveRDS(conserved_markers, "Rdata/seurat_integrated_FindConservedMarkers.rds")