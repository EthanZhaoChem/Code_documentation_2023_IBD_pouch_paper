# R4.1
library(Seurat)
library(ggplot2)
source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')


##################  ################ ##################  ################
# 4. process merged data
seurat <- readRDS('~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged.rds')
path_processed <- '~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged_processed.rds'
plot.dir <- '~/yuzhao1/work/final_RC2rna/0revision/spatial2/plots/'

seurat <- NormalizeData(seurat, assay = 'RNA')
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 400)
seurat <- ScaleData(seurat, vars.to.regress = c('nCount_RNA'))
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunUMAP(seurat,  dims = 1:50, reduction = 'pca', reduction.name = 'umap', 
                  min.dist = 0.2, n.neighbors = 20, seed.use = 1,
                  reduction.key = 'UMAP_')

df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$umap@cell.embeddings)$UMAP_2
df$cluster_name <- seurat$predicted.celltype
p <- plot_df_umap_custom(df, show.label = 'number', custom_colors = custom_colors_rna_all)
png(paste0(plot.dir, 'umap.png'),res = 300, height = 3000, width = 2600)
print(p)
dev.off()

saveRDS(seurat, path_processed)









