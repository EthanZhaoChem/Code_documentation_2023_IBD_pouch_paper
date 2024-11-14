dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')

library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(harmony)
library(Seurat)
library(SeuratObject)
source('~/yuzhao1/scripts/plot.R')
out.dir <- '/project/gca/yuzhao1/work/final_RC2rna/0revision/reviewer3/1control_ethnicity/'


# 1, read controls in pouch paper
seurat_rc2 <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
seurat_rc2 <- subset(seurat_rc2, cells=Cells(seurat_rc2)[seurat_rc2$biopsy_location %in% c('AC', 'TI')])

# 2, read meta data
meta1 <- read.csv('~/yuzhao1/work/final_RC2rna/metadata/update_SampleName.csv', header = T, row.names = NULL)
meta2 <- read.csv('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240311.csv', header = T, row.names = 1)
meta2$Sample_ID <- paste0(meta2$patient, '-', meta2$biopsy_location, '-', meta2$disease_status, meta2$inflammation_status)
meta2$Race[meta2$Race=='Patient Declined'] <- 'Unknown'
meta2 <- meta2[meta2$disease_status== 'Control',]


# # 3, read and process controls in gca paper
# gut_filtered <- readRDS('~/yuzhao1/work/final_GCArna/annotation/rds/gca_combined_final.rds')
# seurat <- subset(gut_filtered, cells=Cells(gut_filtered)[gut_filtered$inflammation_status=='Control'])
# seurat$Sample_ID <- paste0(seurat$Patient_ID,'-',seurat$biopsy_location, '-Control')
# seurat$Race <- mapvalues(seurat$Sample_ID, from = meta2$Sample_ID, to = meta2$Race, warn_missing = F)
# seurat <- NormalizeData(seurat)
# seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
# seurat <- ScaleData(seurat, vars.to.regress = c('percent.mt', 'nCount_RNA', 'CC.Difference'))
# seurat <- RunPCA(seurat, npcs = 50)
# seurat <- RunHarmony(seurat, group.by.vars = 'Patient_ID', max.iter.harmony = 20)
# seurat <- FindNeighbors(seurat, reduction = 'harmony', dims = 1:50)
#  
# seurat <- RunUMAP(seurat,  dims = 1:50, reduction = 'harmony', reduction.name = 'harmony_umap',
#                       min.dist = 0.5, n.neighbors = 50, seed.use = 0,
#                       reduction.key = 'UMAP_')
# saveRDS(seurat, '~/yuzhao1/work/final_RC2rna/0revision/reviewer3/1control_ethnicity/seurat_all_controls.rds')
seurat <- readRDS('~/yuzhao1/work/final_RC2rna/0revision/reviewer3/1control_ethnicity/seurat_all_controls.rds')


# # 4, label transfer
# seurat.anchors <- FindTransferAnchors(reference = seurat_rc2, query = seurat, dims = 1:50, 
#                                       reduction = "pcaproject", verbose = T)
# saveRDS(seurat.anchors, paste0(out.dir, '/seurat.anchors.rds'))

seurat.anchors <- readRDS(paste0(out.dir, '/seurat.anchors.rds'))
predictions <- TransferData(anchorset = seurat.anchors, refdata = seurat_rc2$anno1, dims = 1:50)
seurat <- AddMetaData(seurat, metadata = predictions)

df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$harmony_umap@cell.embeddings)$harmony_umap_1
df$embedding2 <- data.frame(seurat@reductions$harmony_umap@cell.embeddings)$harmony_umap_2
df$cluster_name <- seurat$predicted.id
p<-plot_df_umap_custom(df, show.label = 'number')
png(paste0(out.dir, '/plots/1umap_predicted.png'),res = 300, height = 3000, width = 2600)
print(p)
dev.off()


# 5, labeling
seurat$anno1 <- seurat$predicted.id
seurat$anno1[seurat$anno1 %in% c("cDC1", "cDC2", "Lymphoid DC")] <- 'DC'
seurat$anno1[seurat$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3", "Myofibroblast")] <- 'Fibroblast'
seurat$anno1[seurat$anno1 %in% c("CD4 Tcm" , "Treg", "CD103- CD4 Trm", "CD103+ CD4 Trm")] <- 'CD4T'
seurat$anno1[seurat$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm")] <- 'CD8T'
seurat$anno1[seurat$anno1 %in% c("GC B", "Naive B", "Memory B")] <- 'B'
seurat$anno1[seurat$anno1 %in% c("Monocyte", "Macrophage")] <- 'Macrophage'
seurat$anno1[seurat$anno1 %in% c("Arterial", "Venous", "Lymphatic endothelium")] <- 'Endothelial'
seurat$anno1[seurat$anno1 %in% c("Pericyte", "Contractile pericyte")] <- 'Pericyte'
seurat$anno1[seurat$anno1 %in% c("IgA plasma", "IgG plasma")] <- 'Plasma'
seurat$anno1[seurat$anno1 %in% c('Neutrophil', 'NK T', 'Paneth', 'MAIT', 'ILCs')] <- 'rare'
sort(table(seurat$anno1))
seurat@meta.data[, grepl('prediction', colnames(seurat@meta.data))] <- NULL
seurat@meta.data[, grepl('predicted_anno1', colnames(seurat@meta.data))] <- NULL
seurat@meta.data[, grepl('RNA_snn', colnames(seurat@meta.data))] <- NULL
seurat@meta.data[, grepl('seurat_clusters', colnames(seurat@meta.data))] <- NULL
seurat@meta.data[, grepl('anno2', colnames(seurat@meta.data))] <- NULL

df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$harmony_umap@cell.embeddings)$harmony_umap_1
df$embedding2 <- data.frame(seurat@reductions$harmony_umap@cell.embeddings)$harmony_umap_2
df$cluster_name <- seurat$anno1
p<-plot_df_umap_custom(df, show.label = 'number')
png(paste0(out.dir, '/plots/1umap_anno1.png'),res = 300, height = 3000, width = 2600)
print(p)
dev.off()

saveRDS(seurat, '~/yuzhao1/work/final_RC2rna/0revision/reviewer3/1control_ethnicity/seurat_all_controls.rds')


























