library(seurat)
library(tidyverse)
library(metap)
seurat_integrated_12PC = readRDS("Rdata/seurat_integrated_12PC_res_0.3_Hahn.rds")

#It is tough like that. I could only identify oligodendrocytes. Therefore, need to Identify conserved markers in all conditions. As I have two conditions, need to run FindConservedMarkers
DefaultAssay(seurat_integrated_12PC) <- "RNA"

##Lots of reading to do again in 09_merged_SC_marker_identification.md

# Create function to get conserved markers for any given cluster
annotations = readRDS("Rdata/mouse-gene-annotation.rds")


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
system.time(conserved_markers <- map_dfr(c(0,21), get_conserved))


# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + stim_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

# Visualize top 10 markers per cluster
#View(top10)


saveRDS(conserved_markers, "Rdata/seurat_integrated_12PC_FindConservedMarkers.rds")

