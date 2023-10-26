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
input.dir <- '~/yuzhao1/work/final_RC2atac/peaks/4homer/input/'

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj_tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj_bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj_myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj_others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")

DARs <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_Goblet.rds')


############################ 1. Goblet: ACvsTI ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'Goblet_ACvsTI'
group.1 <- "Goblet-AC"
group.2 <- "Goblet-TI"
annotation_group <- "anno1.loc"
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

df_stat <- archr_helper_markerPeaks_converter(DARs[[contrast_name]])
archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                   start(proj@peakSet),'-',end(proj@peakSet)))
idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
df_stat <- bind_cols(df_stat, peaks_info)


# export DARs in this contrast
df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
df <- data.frame(seqnames = df_stat_significant$seqnames,
                 starts = df_stat_significant$start - 1,
                 ends = df_stat_significant$end, 
                 strands = df_stat_significant$strand)

write.table(df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)


# export DARs in the opposite contrast
df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC < 0-cutoff_Log2FC), ]
df <- data.frame(seqnames = df_stat_significant$seqnames,
                 starts = df_stat_significant$start - 1,
                 ends = df_stat_significant$end, 
                 strands = df_stat_significant$strand)

write.table(df, file=paste0(input.dir, contrast_name, "_opposite.bed"), quote=F, sep="\t", row.names=F, col.names=F)


############################ 2. Goblet: POUvsPP ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'Goblet_POUvsPP'
group.1 <- "Goblet-POU"
group.2 <- "Goblet-PP"
annotation_group <- "anno1.loc"
cutoff_FDR <- 0.01 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 1
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

df_stat <- archr_helper_markerPeaks_converter(DARs[[contrast_name]])
archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                   start(proj@peakSet),'-',end(proj@peakSet)))
idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
df_stat <- bind_cols(df_stat, peaks_info)


# export DARs in this contrast
df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
df <- data.frame(seqnames = df_stat_significant$seqnames,
                 starts = df_stat_significant$start - 1,
                 ends = df_stat_significant$end, 
                 strands = df_stat_significant$strand)

write.table(df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)


# export DARs in the opposite contrast
df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC < 0-cutoff_Log2FC), ]
df <- data.frame(seqnames = df_stat_significant$seqnames,
                 starts = df_stat_significant$start - 1,
                 ends = df_stat_significant$end, 
                 strands = df_stat_significant$strand)

write.table(df, file=paste0(input.dir, contrast_name, "_opposite.bed"), quote=F, sep="\t", row.names=F, col.names=F)



