library(Seurat)
library(arrow)
library(progressr)
library(udunits2)
library(sf)
library(IRIS)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggrepel)
source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')
ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
out.dir <- '~/yuzhao1/work/final_RC2rna/0revision/spatial3/results/'

# 1, prepare xenium input
n_samples <- 4
paths_samples <- c('pouR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__72pov-222__20240711__171029',
                   'ppR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__72pp-222__20240711__171029',
                   'pouR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__72pov-314__20240711__171029',
                   'ppR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__72pp-314__20240711__171029')

xenium.obj.list <- list()
spatial_countMat_list <- list()
spatial_location_list <- list()
for (i in 1:n_samples) {
  sample_name <- names(paths_samples)[[i]]
  tmp.obj <- LoadXenium_alt(paths_samples[[sample_name]], fov = "fov")
  xenium.obj.list[[sample_name]] <- subset(tmp.obj, subset = nCount_Xenium > 0)
  spatial_countMat_list[[sample_name]] <- xenium.obj.list[[sample_name]]@assays$Xenium$counts
  spatial_location_list[[sample_name]] <- xenium.obj.list[[sample_name]]@images$fov$centroids@coords %>% as.data.frame()
  rownames(spatial_location_list[[sample_name]]) <- colnames(spatial_countMat_list[[sample_name]])
}

# 2, prepare single cell input (merge close celltypes)
ref@meta.data$anno1[ref@meta.data$anno1 %in% c('EC1-1', 'EC1-2', 'EC2-1', 'EC2-2')] <- 'EC'
ref@meta.data$anno1[ref@meta.data$anno1 %in% c('Goblet1', 'Goblet2')] <- 'Goblet'
ref@meta.data$anno1[ref@meta.data$anno1 %in% c("CD4 Tcm", "Treg", "CD103- CD4 Trm", "CD103+ CD4 Trm", 
                                               "CD103+ CD8 Trm", "KLRG1+ CD8 Trm", "gdT", "NK T")] <- 'Tcell'
ref@meta.data$anno1[ref@meta.data$anno1 %in% c("Naive B", "GC B", "Memory B")] <- 'Bcell'
ref@meta.data$anno1[ref@meta.data$anno1 %in% c("IgA plasma", "IgG plasma")] <- 'PlasmaB'
ref@meta.data$anno1[ref@meta.data$lineage == 'myeloid'] <- 'myeloid'
ref@meta.data$anno1[ref@meta.data$anno1 %in% c("Myofibroblast", "Stromal-1", "Stromal-3", "Stromal-2")] <- 'Stromal'
ref@meta.data$anno1[ref@meta.data$anno1 %in% c("Arterial", "Pericyte", "Glial", "Venous", "Contractile pericyte",
                                               "Smooth muscle", "Lymphatic endothelium")] <- 'Endothelium'
sc_count <- ref[["RNA"]]$counts
sc_meta <- ref@meta.data
ct.varname <- c("anno1") # Just need to provide the column's name 
sample.varname <- c("Sample_ID") # Provide the column's name for sample sources

# 3, IRIS (run once and save the results)
IRIS_object_raw <- createIRISObject(
  spatial_countMat_list = spatial_countMat_list,
  spatial_location_list = spatial_location_list,
  sc_count = sc_count,
  sc_meta = sc_meta,
  ct.varname = ct.varname,
  sample.varname = sample.varname,
  minCountGene = 20,
  minCountSpot =5) 

IRIS_object_list <- list()
for(numCluster in 3:10){
  clusterID <- paste0('nCluster_', numCluster)
  IRIS_object <- IRIS_spatial(IRIS_object_raw, numCluster = numCluster)
  IRIS_object_list[[clusterID]] <- IRIS_object
}
saveRDS(IRIS_object_list, paste0(out.dir, 'IRIS_object_list.rds'))

# 4, plot domain location
IRIS_object_list <- readRDS(paste0(out.dir, 'IRIS_object_list.rds'))
colors <- c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','black')
for(numCluster in 3:10){
  clusterID <- paste0('nCluster_', numCluster)
  IRIS_object <- IRIS_object_list[[clusterID]] 
  domain_results <-  IRIS_object@spatialDomain[,c("Slice","spotName","IRIS_domain")]
  spatial_location <- IRIS_object@spatialDomain[,c("x","y")]
  p1 <- IRIS.visualize.domain(domain_results, spatial_location, colors = colors, numCols = 4)
  pdf(paste0(out.dir, 'domain_location/spatial_domain_', numCluster, '.pdf'), width = 20, height = 6)
  print(p1)
  dev.off()
}

# 5, differential gene for different domains, heatmap
seurat_xe <- readRDS('~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged.rds')
seurat_xe <- ScaleData(seurat_xe, vars.to.regress = c('nCount_RNA'), features = rownames(seurat)) # to help plot
seurat_xe$cellID <- Cells(seurat_xe)

for(numCluster in 3:10){
  clusterID <- paste0('nCluster_', numCluster)
  IRIS_object <- IRIS_object_list[[clusterID]] 
  domain_results <-  IRIS_object@spatialDomain[,c("Slice","spotName","IRIS_domain")]
  domain_results$cellID <- paste0(domain_results$Slice, '_', domain_results$spotName)
  
  seurat <- seurat_xe
  seurat <- subset(seurat, cells = Cells(seurat)[seurat$cellID %in% domain_results$cellID])
  seurat$domain <- mapvalues(seurat$cellID, domain_results$cellID, domain_results$IRIS_domain, warn_missing = F)
  Idents(seurat) <- seurat$domain
  all.markers <- FindAllMarkers(seurat, only.pos = TRUE)
  
  all.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  
  pdf(paste0(out.dir, 'domain_deg_heatmap/spatial_domain_', numCluster, '_deg_heatmap.pdf'), width = 6, height = numCluster * 2)
  p2 <- DoHeatmap(seurat, features = top10$gene) + NoLegend()
  print(p2)
  dev.off()
}














