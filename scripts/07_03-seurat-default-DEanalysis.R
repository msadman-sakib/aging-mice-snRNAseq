library(Seurat)
library(cowplot)
library(patchwork)
library(tidyverse)

# Bring in Seurat object
seurat = readRDS("Rdata/seurat_integrated_cell_labelled_final.rds")

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





