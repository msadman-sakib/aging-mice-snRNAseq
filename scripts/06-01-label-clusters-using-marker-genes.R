#Now, loading the outputs from HPC outputs
library(Seurat)
library(tidyverse)
library(rio)
library(enrichR)
library(biomaRt)

#load seurat integrated dataset.
seurat_integrated = readRDS("Rdata/seurat_integrated_PC40_res0.3.rds")

# Like in my 05_clustering-QC.R script, I am gonna do the marker gene checking here. but not for individual genes, but for genelist (with AddModuleScore() function)


filelist = list.files("marker_genelist", full.names = T)
marker.genelist = lapply(filelist, function(x)read_table(x,col_names = F))
names(marker.genelist) = tools::file_path_sans_ext(basename(filelist))
for(i in 1:length(marker.genelist)){
  marker.genelist[[i]] = marker.genelist[[i]] %>%  pull()
}

#Now biomartr
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")#, host="useast.ensembl.org") ##changing the host to usa made it work! How weird!

# #See the filters to choose from. This is input data, can be for example, ensemble gene id
# filters = listFilters(ensembl)
# filters[1:5,]
# 
# #This is output data. take a look to find
# attributes = listAttributes(ensembl)
# attributes[1:5,]

for(i in 1:length(marker.genelist)){
marker.genelist[[i]] = getBM(filters = 'ensembl_gene_id',
                              attributes=c('external_gene_name'), #multiple attributes can be selected here. 
                              values = marker.genelist[[i]],
                              mart = ensembl) %>% pull(external_gene_name)
}


##For AddModuleScore fuction
DefaultAssay(seurat_integrated) = "RNA"

seurat_integrated <- AddModuleScore(object = seurat_integrated, features = marker.genelist[1], ctrl = 5, name = names(marker.genelist)[1])
names(x = seurat_integrated[[]])
FeaturePlot(object = seurat_integrated, features = paste0(names(marker.genelist)[1],"1"), reduction = "umap", order = TRUE, min.cutoff = 'q10', label = TRUE)





# # Rename all identities
# seurat_integrated <- RenameIdents(object = seurat_integrated, 
#                                   "0"="0 Glutamatergic DG 1",
#                                   "1"="1 Glutamatergic CA1-do",
#                                   "2"="2 Oligo",
#                                   "3"="3 GABAergic Sncg 1",
#                                   "4"="4 GABAergic Ntng1 HPF 1",
#                                   "5"="5 Glutamatergic SUB",
#                                   "6"="6 GABAergic Pax6",
#                                   "7"="7 Glutamatergic CA3-do",
#                                   "8"="8 Oligo-OPC 1",
#                                   "9"="9 GABAergic Sncg 2",
#                                   "10"="10 Glutamatergic DG 2",
#                                   "11"="11 Glutamatergic NP SUB",
#                                   "12"="12 Astro",
#                                   "13"="13 Glutamatergic DG 3",
#                                   "14"="14 Glutamatergic L2 IT ENTl",
#                                   "15"="15 GABAergic Ntng1 HPF 2",
#                                   "16"="16 Glutamatergic Car3",
#                                   "17"="17 Oligo-Glutamatergic",
#                                   "18"="18 Micro 1",
#                                   "19"="19 GABAergic Lamp5",
#                                   "20"="20 GABAergic Sst",
#                                   "21"="21 Glutamatergic Mossy",
#                                   "22"="22 Glutamatergic L2 IT APr",
#                                   "23"="23 Glutamatergic CR",
#                                   "24"="24 Gabaergic Glutamatergic",
#                                   "25"="25 GABAergic Sst Lamp5 Lhx6",
#                                   "26"="26 Peri Endo",
#                                   "27"="27 Micro 2",
#                                   "28"="28 Micro 3",
#                                   "29"="29 Micro 4",
#                                   "30"="30 Oligo 3",
#                                   "31"="31 Oligo-OPC 2")
# 
# # Plot the UMAP
# DimPlot(object = seurat_integrated, 
#         reduction = "umap", 
#         label = TRUE,
#         label.size = 2,
#         repel = TRUE)
# ggsave("plots/final-plots/UMAP_clusters_relabelled.pdf", height = 6, width = 12)
# 
# saveRDS(seurat_integrated, "Rdata/seurat_integrated_cell_labelled.rds")
# 
# 
