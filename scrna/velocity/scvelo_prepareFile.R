dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')

# customize 
seurat_obj <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/epithelial3.rds')
out_data_dir <- '~/yuzhao1/work/final_RC2rna/velocity/epithelial/'

# output name
out_metadata_file <- paste0(out_data_dir, 'metadata.csv')
out_ExprMtx_file <- paste0(out_data_dir, 'counts.mtx')
out_DimRec_file <- paste0(out_data_dir, 'DimRec.csv')
out_genes_file <- paste0(out_data_dir, 'gene_names.csv')

# write meta.data
df <- data.frame(seurat_obj@meta.data)
df$embedding1 <- data.frame(seurat_obj@reductions$harmony_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat_obj@reductions$harmony_umap@cell.embeddings)$UMAP_2
df$barcode <- rownames(seurat_obj@meta.data)
write.csv(df, file=out_metadata_file, quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=out_ExprMtx_file)

# write dimesnionality reduction matrix, in this example case harmony matrix
write.csv(seurat_obj@reductions$harmony@cell.embeddings, file=out_DimRec_file, quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file=out_genes_file ,
  quote=F,row.names=F,col.names=F
)

seurat_obj$anno1[which(seurat_obj$anno1 == "EC1-MMP1")] <- 'EC1-3'

DimPlot(seurat_obj, reduction = 'harmony_umap', group.by = 'anno1')


