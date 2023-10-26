dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/helper_archr.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')
addArchRThreads(4)
tf_names <- readRDS('~/yuzhao1/work/manu/rc2/plots/6tf_logo/tf_names.rds')

# 0. build the rna filtering matrix (tf by anno2)
seurat <- readRDS("~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds")
seurat.list <- list()
for (temp.anno in unique(seurat$anno2)){
  seurat.list[[temp.anno]] <- subset(seurat, anno2 == temp.anno)
}

genes.broadlyExpressed.pool <- list()
for (temp.anno in unique(seurat$anno2)){
  seurat.temp <- seurat.list[[temp.anno]]
  ncellsPerGene <- rowSums2(seurat.temp@assays$RNA@counts > 0.5) # cutoff is 0.5 gene for a cell (consider >0.5 as 1)
  genes.broadlyExpressed.idx <- order(ncellsPerGene/ncol(seurat.temp), decreasing = T)[1:3000] # top5000
  temp <- rownames(seurat.temp)[genes.broadlyExpressed.idx]
  genes.broadlyExpressed.pool[[temp.anno]] <- temp
}



# 1. epithelial
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
DARs <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno2_epithelial', '.rds'))
DARs_plot <- list()
df_customized_enrichment <- tf_names
df_customized_enrichment.list <- list()
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1

df_stat_all <- archr_helper_markerPeaks_converter_multiple(DARs)
for(cluster.name in names(df_stat_all)){
  df_stat <- df_stat_all[[cluster.name]]
  archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                     start(proj@peakSet),'-',end(proj@peakSet)))
  idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
  peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
  df_stat <- bind_cols(df_stat, peaks_info)

  # significant peaks
  df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
  DARs_plot[[cluster.name]] <- rownames(df_stat_significant)

  # enrich for motifs
  usefulpeaks_vector <- DARs_plot[[cluster.name]]
  df_customized_enrichment.list[[cluster.name]] <- archr_customized_motif_enrichment(proj,
                                                                                 peakAnnotation = "Motif",
                                                                                 candidate_peaks_vector = usefulpeaks_vector)

  # save the mlog10Padj column based on the same motif order
  df_customized_enrichment[[cluster.name]] <- df_customized_enrichment.list[[cluster.name]][df_customized_enrichment$tf_full, 'mlog10Padj']
}

saveRDS(df_customized_enrichment, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno2_epithelial_enriched_mlog10Padj', '.rds'))

df_customized_enrichment <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno2_epithelial_enriched_mlog10Padj', '.rds'))
setdiff(unique(proj$anno2), unique(seurat$anno2)) # make sure each atac celltype is represented

tf_filtered.list <- list()
for (celltype in unique(proj$anno2)) {
  temp_celltype <- df_customized_enrichment[, c('tf', celltype)]
  temp_celltype <- temp_celltype[temp_celltype[[celltype]] > 3,]
  tf.top20 <- temp_celltype$tf[order(temp_celltype[[celltype]], decreasing = T)] %>% head(20)
  tf_filtered.list[[celltype]] <- tf.top20[tf.top20 %in% genes.broadlyExpressed.pool[[celltype]]]
  df_customized_enrichment
}
tf_filtered <- unique(unlist(tf_filtered.list ))

# set the non-significant values to zero (induced by TFs in other celltypes)
df_filtered <- df_customized_enrichment[df_customized_enrichment$tf%in%tf_filtered,]
for (i in 1:nrow(df_filtered)) {
  for(j in 3:ncol(df_filtered)){
    df_filtered[i,j] <- ifelse(df_filtered[i,j]<3, 0, df_filtered[i,j])
  }
}
saveRDS(df_filtered ,
        paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno2_epithelial_enriched_mlog10Padj_filtered', '.rds'))





# 2. tcell
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
DARs <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))
DARs_plot <- list()
df_customized_enrichment <- tf_names
df_customized_enrichment.list <- list()
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1

df_stat_all <- archr_helper_markerPeaks_converter_multiple(DARs)
for(cluster.name in names(df_stat_all)){
  df_stat <- df_stat_all[[cluster.name]]
  archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                     start(proj@peakSet),'-',end(proj@peakSet)))
  idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
  peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
  df_stat <- bind_cols(df_stat, peaks_info)

  # significant peaks
  df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
  DARs_plot[[cluster.name]] <- rownames(df_stat_significant)

  # enrich for motifs
  usefulpeaks_vector <- DARs_plot[[cluster.name]]
  df_customized_enrichment.list[[cluster.name]] <- archr_customized_motif_enrichment(proj,
                                                                                 peakAnnotation = "Motif",
                                                                                 candidate_peaks_vector = usefulpeaks_vector)

  # save the mlog10Padj column based on the same motif order
  df_customized_enrichment[[cluster.name]] <- df_customized_enrichment.list[[cluster.name]][df_customized_enrichment$tf_full, 'mlog10Padj']
}

saveRDS(df_customized_enrichment, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))

df_customized_enrichment <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))
setdiff(unique(proj$anno1), unique(seurat$anno1)) # make sure each atac celltype is represented

tf_filtered.list <- list()
for (celltype in unique(proj$anno1)) {
  temp_celltype <- df_customized_enrichment[, c('tf', celltype)]
  temp_celltype <- temp_celltype[temp_celltype[[celltype]] > 3,]
  tf.top20 <- temp_celltype$tf[order(temp_celltype[[celltype]], decreasing = T)] %>% head(20)
  tf_filtered.list[[celltype]] <- tf.top20[tf.top20 %in% genes.broadlyExpressed.pool[[celltype]]]
  df_customized_enrichment
}
tf_filtered <- unique(unlist(tf_filtered.list))

# set the non-significant values to zero (induced by TFs in other celltypes)
df_filtered <- df_customized_enrichment[df_customized_enrichment$tf%in%tf_filtered,]
for (i in 1:nrow(df_filtered)) {
  for(j in 3:ncol(df_filtered)){
    df_filtered[i,j] <- ifelse(df_filtered[i,j]<3, 0, df_filtered[i,j])
  }
}

saveRDS(df_filtered,
        paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj_filtered', '.rds'))






# 3. bcell
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
DARs <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))
DARs_plot <- list()
df_customized_enrichment <- tf_names
df_customized_enrichment.list <- list()
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1

df_stat_all <- archr_helper_markerPeaks_converter_multiple(DARs)
for(cluster.name in names(df_stat_all)){
  df_stat <- df_stat_all[[cluster.name]]
  archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                     start(proj@peakSet),'-',end(proj@peakSet)))
  idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
  peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
  df_stat <- bind_cols(df_stat, peaks_info)
  
  # significant peaks
  df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
  DARs_plot[[cluster.name]] <- rownames(df_stat_significant)
  
  # enrich for motifs
  usefulpeaks_vector <- DARs_plot[[cluster.name]]
  df_customized_enrichment.list[[cluster.name]] <- archr_customized_motif_enrichment(proj,
                                                                                     peakAnnotation = "Motif",
                                                                                     candidate_peaks_vector = usefulpeaks_vector)
  
  # save the mlog10Padj column based on the same motif order
  df_customized_enrichment[[cluster.name]] <- df_customized_enrichment.list[[cluster.name]][df_customized_enrichment$tf_full, 'mlog10Padj']
}

saveRDS(df_customized_enrichment, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))

df_customized_enrichment <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))

genes.broadlyExpressed.pool[["Plasma cell"]] <- unique(c(genes.broadlyExpressed.pool$`IgA plasma`, genes.broadlyExpressed.pool$`IgG plasma`))
setdiff(unique(proj$anno1), names(genes.broadlyExpressed.pool)) # make sure each atac celltype is represented


tf_filtered.list <- list()
for (celltype in unique(proj$anno1)) {
  temp_celltype <- df_customized_enrichment[, c('tf', celltype)]
  temp_celltype <- temp_celltype[temp_celltype[[celltype]] > 3,]
  tf.top20 <- temp_celltype$tf[order(temp_celltype[[celltype]], decreasing = T)] %>% head(20)
  tf_filtered.list[[celltype]] <- tf.top20[tf.top20 %in% genes.broadlyExpressed.pool[[celltype]]]
  df_customized_enrichment
}
tf_filtered <- unique(unlist(tf_filtered.list))

# set the non-significant values to zero (induced by TFs in other celltypes)
df_filtered <- df_customized_enrichment[df_customized_enrichment$tf%in%tf_filtered,]
for (i in 1:nrow(df_filtered)) {
  for(j in 3:ncol(df_filtered)){
    df_filtered[i,j] <- ifelse(df_filtered[i,j]<3, 0, df_filtered[i,j])
  }
}

saveRDS(df_filtered,
        paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj_filtered', '.rds'))





# 4. myeloid
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
DARs <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))
DARs_plot <- list()
df_customized_enrichment <- tf_names
df_customized_enrichment.list <- list()
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1

df_stat_all <- archr_helper_markerPeaks_converter_multiple(DARs)
for(cluster.name in names(df_stat_all)){
  df_stat <- df_stat_all[[cluster.name]]
  archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                     start(proj@peakSet),'-',end(proj@peakSet)))
  idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
  peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
  df_stat <- bind_cols(df_stat, peaks_info)
  
  # significant peaks
  df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
  DARs_plot[[cluster.name]] <- rownames(df_stat_significant)
  
  # enrich for motifs
  usefulpeaks_vector <- DARs_plot[[cluster.name]]
  df_customized_enrichment.list[[cluster.name]] <- archr_customized_motif_enrichment(proj,
                                                                                     peakAnnotation = "Motif",
                                                                                     candidate_peaks_vector = usefulpeaks_vector)
  
  # save the mlog10Padj column based on the same motif order
  df_customized_enrichment[[cluster.name]] <- df_customized_enrichment.list[[cluster.name]][df_customized_enrichment$tf_full, 'mlog10Padj']
}

saveRDS(df_customized_enrichment, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))

df_customized_enrichment <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))

genes.broadlyExpressed.pool[["DC"]] <- unique(c(genes.broadlyExpressed.pool[['cDC1']], 
                                                genes.broadlyExpressed.pool[['cDC2']],
                                                genes.broadlyExpressed.pool[["Lymphoid DC"]]))
setdiff(unique(proj$anno1), names(genes.broadlyExpressed.pool)) # make sure each atac celltype is represented


tf_filtered.list <- list()
for (celltype in unique(proj$anno1)) {
  temp_celltype <- df_customized_enrichment[, c('tf', celltype)]
  temp_celltype <- temp_celltype[temp_celltype[[celltype]] > 3,]
  tf.top20 <- temp_celltype$tf[order(temp_celltype[[celltype]], decreasing = T)] %>% head(20)
  tf_filtered.list[[celltype]] <- tf.top20[tf.top20 %in% genes.broadlyExpressed.pool[[celltype]]]
  df_customized_enrichment
}
tf_filtered <- unique(unlist(tf_filtered.list))

# set the non-significant values to zero (induced by TFs in other celltypes)
df_filtered <- df_customized_enrichment[df_customized_enrichment$tf%in%tf_filtered,]
for (i in 1:nrow(df_filtered)) {
  for(j in 3:ncol(df_filtered)){
    df_filtered[i,j] <- ifelse(df_filtered[i,j]<3, 0, df_filtered[i,j])
  }
}

saveRDS(df_filtered,
        paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj_filtered', '.rds'))





# 5, stromal
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")
DARs <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))
DARs_plot <- list()
df_customized_enrichment <- tf_names
df_customized_enrichment.list <- list()
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1

df_stat_all <- archr_helper_markerPeaks_converter_multiple(DARs)
for(cluster.name in names(df_stat_all)){
  df_stat <- df_stat_all[[cluster.name]]
  archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                     start(proj@peakSet),'-',end(proj@peakSet)))
  idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
  peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
  df_stat <- bind_cols(df_stat, peaks_info)
  
  # significant peaks
  df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
  DARs_plot[[cluster.name]] <- rownames(df_stat_significant)
  
  # enrich for motifs
  usefulpeaks_vector <- DARs_plot[[cluster.name]]
  df_customized_enrichment.list[[cluster.name]] <- archr_customized_motif_enrichment(proj,
                                                                                     peakAnnotation = "Motif",
                                                                                     candidate_peaks_vector = usefulpeaks_vector)
  
  # save the mlog10Padj column based on the same motif order
  df_customized_enrichment[[cluster.name]] <- df_customized_enrichment.list[[cluster.name]][df_customized_enrichment$tf_full, 'mlog10Padj']
}

saveRDS(df_customized_enrichment, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))

df_customized_enrichment <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj', '.rds'))

setdiff(unique(proj$anno1), names(genes.broadlyExpressed.pool)) # make sure each atac celltype is represented

tf_filtered.list <- list()
for (celltype in unique(proj$anno1)) {
  temp_celltype <- df_customized_enrichment[, c('tf', celltype)]
  temp_celltype <- temp_celltype[temp_celltype[[celltype]] > 3,]
  tf.top20 <- temp_celltype$tf[order(temp_celltype[[celltype]], decreasing = T)] %>% head(20)
  tf_filtered.list[[celltype]] <- tf.top20[tf.top20 %in% genes.broadlyExpressed.pool[[celltype]]]
  df_customized_enrichment
}
tf_filtered <- unique(unlist(tf_filtered.list))

# set the non-significant values to zero (induced by TFs in other celltypes)
df_filtered <- df_customized_enrichment[df_customized_enrichment$tf%in%tf_filtered,]
for (i in 1:nrow(df_filtered)) {
  for(j in 3:ncol(df_filtered)){
    df_filtered[i,j] <- ifelse(df_filtered[i,j]<3, 0, df_filtered[i,j])
  }
}

saveRDS(df_filtered,
        paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_enriched_mlog10Padj_filtered', '.rds'))









