library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')
source('~/yuzhao1/scripts/plot.R')

coords.dir <- '~/yuzhao1/work/final_RC2rna/0revision/spatial2/coords/'
seurat <- readRDS('~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged_processed.rds')
seurat$cluster <- seurat$predicted.celltype
metadata <- seurat@meta.data

# 1. cell type label
sample.names <- unique(seurat$sample)
n_sample <- length(sample.names)

for (i in 1:n_sample) {
  tmp.sample <- sample.names[[i]]
  tmp.meta <- metadata[metadata$sample == tmp.sample, ]
  cellnames <- rownames(tmp.meta) %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,2)
  coords <- data.frame(cell_id = cellnames, group = tmp.meta$cluster)
  write.csv(coords, paste0(coords.dir, tmp.sample, '.csv'), row.names = F, quote = F)
}



# 2. gene module