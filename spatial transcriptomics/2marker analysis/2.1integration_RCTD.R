# R4.3
library(Seurat)
library(arrow)
library(progressr)
library(udunits2)
library(sf)
library(ggplot2)
library(spacexr)
source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')

n_samples <- 6
paths_samples <- c('ac' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__50AC__20240711__171029',
                   'ti' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__50TI__20240711__171029',
                   'pouR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__72pov-222__20240711__171029',
                   'ppR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__72pp-222__20240711__171029',
                   'pouR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__72pov-314__20240711__171029',
                   'ppR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__72pp-314__20240711__171029')
# sample 53 can't be mapped because of extremely low quality

xenium.obj.list <- list()
for (i in 1:n_samples) {
  sample_name <- names(paths_samples)[[i]]
  tmp.obj <- LoadXenium_alt(paths_samples[[sample_name]], fov = "fov")
  xenium.obj.list[[sample_name]] <- subset(tmp.obj, subset = nCount_Xenium > 0)
}

##################  ################ ##################  ################
# 1. build reference data
seurat_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
counts <- GetAssayData(seurat_ref, assay = "RNA", slot = "counts")
cluster <- as.factor(seurat_ref$anno1)
names(cluster) <- colnames(seurat_ref)
nUMI <- seurat_ref$nCount_RNA
names(nUMI) <- colnames(seurat_ref)
nUMI <- colSums(counts)
reference <- Reference(counts, cluster, nUMI, require_int = F, n_max_cells = 500, min_UMI = 10)

# 2. integration for each sample
for (i in 1:n_samples) {
  # calculation
  sample_name <- names(paths_samples)[[i]]
  xenium.obj <- xenium.obj.list[[sample_name]]
  query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")
  coords <- GetTissueCoordinates(xenium.obj, which = "centroids")
  rownames(coords) <- coords$cell
  coords$cell <- NULL
  query <- SpatialRNA(coords, query.counts, colSums(query.counts))
  RCTD <- create.RCTD(query, reference, max_cores = 10, DOUBLET_THRESHOLD = 20, 
                      counts_MIN=5, CELL_MIN_INSTANCE = 20, UMI_min = 10)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  # labeling and subseting
  annotations.df <- RCTD@results$results_df
  annotations <- annotations.df$first_type
  names(annotations) <- rownames(annotations.df)
  xenium.obj$predicted.celltype <- annotations
  keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
  xenium.obj.list[[sample_name]] <- subset(xenium.obj, cells = keep.cells)
}

saveRDS(xenium.obj.list, '~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.obj.list.rds')



##################  ################ ##################  ################
# 3. merge data
xenium.obj <- merge(xenium.obj.list[[1]], 
                    y = xenium.obj.list[2:n_samples], 
                    add.cell.ids = names(paths_samples), 
                    project = "rc2")

xenium.obj$sample <- sapply(strsplit(colnames(xenium.obj), "_"), `[`, 1)
xenium.obj[['Xenium']] <- JoinLayers(xenium.obj[['Xenium']])

xenium.merged <- CreateSeuratObject(counts = xenium.obj@assays$Xenium$counts,
                                    meta.data = xenium.obj@meta.data)

# make it compatible for R4.0
xenium.merged[["RNA"]] <- as(object = xenium.merged[["RNA"]], Class = "Assay")
saveRDS(xenium.merged, '~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged.rds')





