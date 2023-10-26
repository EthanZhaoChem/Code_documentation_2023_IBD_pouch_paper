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
# our.dir <- '~/yuzhao1/work/manu/rc2/plots/4atac_MarkerPeaks/'

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj_tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj_bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj_myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj_others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")

DARs <- list()
DARs_significant <- list()
motif_stats <- list()


############################ 1. Goblet: AC vs TI ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'Goblet_ACvsTI'
group.1 <- "Goblet-AC"
group.2 <- "Goblet-TI"
annotation_group <- "anno1.loc"
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1

# calculation
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
DARs[[contrast_name]] <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = annotation_group,
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = group.1,
  bgdGroups = group.2
)

# clean up data
df_stat <- archr_helper_markerPeaks_converter(DARs[[contrast_name]])
archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                   start(proj@peakSet),'-',end(proj@peakSet)))
idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
df_stat <- bind_cols(df_stat, peaks_info)

df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & abs(df_stat$Log2FC) > cutoff_Log2FC), ]
DARs_significant[[contrast_name]] <- df_stat_significant

# motifs up
motifsUp <- peakAnnoEnrichment(
  seMarker = DARs[[contrast_name]],
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC > 0.5"
)
motif_stat_up <- archr_helper_markerTFs_converter(motifsUp)
motif_stats[[paste0(contrast_name, '_up')]] <- motif_stat_up 

# motifs down
motifsDown <- peakAnnoEnrichment(
  seMarker = DARs[[contrast_name]],
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC < -0.5"
)
motif_stat_down <- archr_helper_markerTFs_converter(motifsDown)
motif_stats[[paste0(contrast_name, '_down')]] <- motif_stat_down




############################ 2. Goblet: POU vs PP ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'Goblet_POUvsPP'
group.1 <- "Goblet-POU"
group.2 <- "Goblet-PP"
annotation_group <- "anno1.loc"
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1

# calculation
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
DARs[[contrast_name]] <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = annotation_group,
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = group.1,
  bgdGroups = group.2
)

# clean up data
df_stat <- archr_helper_markerPeaks_converter(DARs[[contrast_name]])
archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                   start(proj@peakSet),'-',end(proj@peakSet)))
idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
df_stat <- bind_cols(df_stat, peaks_info)

df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & abs(df_stat$Log2FC) > cutoff_Log2FC), ]
DARs_significant[[contrast_name]] <- df_stat_significant

# motifs up
motifsUp <- peakAnnoEnrichment(
  seMarker = DARs[[contrast_name]],
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC > 0.5"
)
motif_stat_up <- archr_helper_markerTFs_converter(motifsUp)
motif_stats[[paste0(contrast_name, '_up')]] <- motif_stat_up 

# motifs down
motifsDown <- peakAnnoEnrichment(
  seMarker = DARs[[contrast_name]],
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC < -0.5"
)
motif_stat_down <- archr_helper_markerTFs_converter(motifsDown)
motif_stats[[paste0(contrast_name, '_down')]] <- motif_stat_down




# ############################ 0, summary ###############################
saveRDS(DARs, '~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_Goblet.rds')
saveRDS(DARs_significant, '~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_significant_Goblet.rds')
saveRDS(motif_stats, '~/yuzhao1/work/final_RC2atac/peaks/3DARs/motif_stats_Goblet.rds')

# DARs <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_Goblet.rds')
# DARs_significant <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_significant_Goblet.rds')
# motif_stats <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/3DARs/motif_stats_Goblet.rds')
# 
# contrast_name <- 'EC_ACvsTI'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# rna_deg_markers_ec_tfs %>% 
#   intersect(., table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100) %>% names())
# 
# contrast_name <- 'EC_POU2vsPOU1'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# rna_deg_markers_ec_tfs %>% 
#   intersect(., table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100) %>% names())
# 
# contrast_name <- 'EC_POU2vsPP'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# 
# contrast_name <- 'EC_POU1vsPP'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# 
# contrast_name <- 'EC_POUvsPP'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# 
# contrast_name <- 'Stem_ACvsTI'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# 
# contrast_name <- 'Stem_POUvsPP'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# 
# contrast_name <- 'Goblet_ACvsTI'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)
# 
# contrast_name <- 'Goblet_POUvsPP'
# table(DARs_significant[[contrast_name]]$nearestGene) %>% sort(decreasing = T) %>% head(100)





