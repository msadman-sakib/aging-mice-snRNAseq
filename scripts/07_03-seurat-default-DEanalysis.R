library(Seurat)
library(cowplot)
library(patchwork)
library(tidyverse)
library(DESeq2)
library(MAST)

# Bring in Seurat object
seurat = readRDS("Rdata/seurat_integrated_cell_labelled_final.rds")
DefaultAssay(seurat) = "RNA"
# let's check average expression between conditions in all cell types. 
theme_set(theme_cowplot())

##It can be done in a for loop.
# for(i in levels(Idents(seurat))){
# subset.cell <- subset(seurat, idents = i)
# Idents(subset.cell) <- "sample"
# avg.subset.cell <- log1p(AverageExpression(subset.cell, verbose = T)$RNA)
# ggplot(as.data.frame(avg.subset.cell), aes(young, old)) + geom_point() + ggtitle(i)
# ggsave(paste0("DE_analysis_seurat_default/figures/",i,".pdf"),height = 5, width = 5)
# }
# 

#But, lets write a function and use map() from purrr to do this. This is wayyyy faster than a for loop.
avg.expression.plot.fn = function(x){subset.cell <- subset(seurat, idents = x)
                                  Idents(subset.cell) <- "sample"
                                  avg.subset.cell <- log1p(AverageExpression(subset.cell, verbose = T)$RNA)
                                  ggplot(as.data.frame(avg.subset.cell), aes(young, old)) + geom_point() + ggtitle(x)
                                  ggsave(paste0("DE_analysis_seurat_default/figures/",x,".pdf"),height = 5, width = 5)
                                  }
map(levels(Idents(seurat)),avg.expression.plot.fn)

#----------------Differential expression using Seurat FindMarkers----------------#
seurat.copy = seurat #making a copy not to manipulate the original seurat object
DefaultAssay(seurat.copy) = "RNA"
seurat.copy$celltype.sample <- paste(Idents(seurat.copy), seurat.copy$sample, sep="_")
seurat.copy$celltype <- Idents(seurat.copy)
Idents(seurat.copy) <- "celltype.sample"

##Running default differential expression (using wilcoxon)
seurat.diff.expression.default.fn = function(x){ident1 <- paste0(x,"_old")
                                                ident2 <- paste0(x,"_young")
                                                sample.diffgenes <- FindMarkers(seurat.copy, ident.1 = ident1, ident.2=ident2)
                                                write.csv(sample.diffgenes, file=paste0("DE_analysis_seurat_default/results/",x,"_oldvsyoung",".csv"))
                                                }
map(levels(seurat.copy$celltype),seurat.diff.expression.default.fn)

#NOTE: I observed huge numbers of diff. expressed genes with super significance. That was expected like I saw before. Now, I am gonna try two other algorithms. First, MAST and then DESeq2. MAST might be the way to go, but lets see. 

#DE testing using MAST
seurat.diff.expression.MAST.fn = function(x){ident1 <- paste0(x,"_old")
                                                ident2 <- paste0(x,"_young")
                                                sample.diffgenes <- FindMarkers(seurat.copy, ident.1 = ident1, ident.2=ident2, test.use = "MAST")
                                                write.csv(sample.diffgenes, file=paste0("DE_analysis_seurat_default/results-mast-test/",x,"_oldvsyoung",".csv"))
                                                }
map(levels(seurat.copy$celltype),seurat.diff.expression.MAST.fn)

#DE testing using DESeq2, but takes a lot time. Run this chunk in HPC.
seurat.diff.expression.DESEQ2.fn = function(x){ident1 <- paste0(x,"_old")
                                                ident2 <- paste0(x,"_young")
                                                sample.diffgenes <- FindMarkers(seurat.copy, ident.1 = ident1, ident.2=ident2, test.use = "DESeq2")
                                                write.csv(sample.diffgenes, file=paste0("DE_analysis_seurat_default/results-deseq2-test/",x,"_oldvsyoung",".csv"))
                                                }
map(levels(seurat.copy$celltype),seurat.diff.expression.DESEQ2.fn)



############# to be run in HPC #########
# library(Seurat)
# library(DESeq2)
# library(tidyverse)
# 
# seurat = readRDS("seurat_integrated_cell_labelled_final.rds")
# DefaultAssay(seurat) = "RNA"
# seurat$celltype.sample <- paste(Idents(seurat), seurat$sample, sep="_")
# seurat$celltype <- Idents(seurat)
# Idents(seurat) <- "celltype.sample"
# 
# seurat.diff.expression.DESEQ2.fn = function(x){ident1 <- paste0(x,"_old")
# ident2 <- paste0(x,"_young")
# sample.diffgenes <- FindMarkers(seurat, ident.1 = ident1, ident.2=ident2, test.use = "DESeq2")
# write.csv(sample.diffgenes, file=paste0("results/",x,"_oldvsyoung",".csv"))
# }
# map(levels(seurat$celltype),seurat.diff.expression.DESEQ2.fn)




