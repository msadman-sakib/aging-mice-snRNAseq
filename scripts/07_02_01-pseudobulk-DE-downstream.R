library(tidyverse)


tbl.names <- list.files(path = "DE_analysis_scrnaseq/results/",pattern = "*all_genes.csv")
tbl.names = tools::file_path_sans_ext(tbl.names)


tbl <- list.files(path = "DE_analysis_scrnaseq/results/",pattern = "*all_genes.csv", full.names = T) %>% map(read_csv)

names(tbl) = tbl.names


tbl.downregulated = lapply(tbl,function(x){x %>% filter(padj < 0.1, log2FoldChange < -0.26)})
tbl.downregulated = data.frame(lapply(tbl.downregulated,function(x){x %>% nrow()}))
tbl.downregulated$index = c(1)
tbl.downregulated = pivot_longer(tbl.downregulated,!index, names_to = "cell.types", values_to = "upregulated")
tbl.downregulated = tbl.downregulated %>% select(-index)

tbl.upregulated = lapply(tbl,function(x){x %>% filter(padj < 0.1, log2FoldChange > 0.26)})
tbl.upregulated = data.frame(lapply(tbl.upregulated,function(x){x %>% nrow()}))
tbl.upregulated$index = c(1)
tbl.upregulated = pivot_longer(tbl.upregulated,!index, names_to = "cell.types", values_to = "downregulated")
tbl.upregulated = tbl.upregulated %>% select(-index)

tbl.up.down.summary = inner_join(tbl.upregulated,tbl.downregulated, by = "cell.types")
tbl.up.down.summary = pivot_longer(tbl.up.down.summary , !cell.types, names_to = "changes", values_to = "count") 

tbl.up.down.summary.notzeros = tbl.up.down.summary %>% filter(!count==0)

tbl.up.down.summary.notzeros$cell.types = str_remove(tbl.up.down.summary.notzeros$cell.types, "_old_vs_young_all_genes")

#Stacked bar plot showing the number of DE genes detected in each cell type
ggplot(tbl.up.down.summary.notzeros, aes(fill=cell.types, y=count, x=changes)) + 
  geom_bar(position="stack", stat="identity") + xlab("") + scale_fill_brewer(palette="Paired")
