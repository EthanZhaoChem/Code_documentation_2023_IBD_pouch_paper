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



DARs <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_Stem.rds')


############################ 1. Stem: ACvsTI ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'Stem_ACvsTI'
group.1 <- "Stem-AC"
group.2 <- "Stem-TI"
annotation_group <- "anno1.loc"
cutoff_FDR <- 0.1 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 0.5
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

# export DARs in this contrast
df_stat <- archr_helper_markerPeaks_converter(DARs[[contrast_name]])
archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                   start(proj@peakSet),'-',end(proj@peakSet)))
idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
df_stat <- bind_cols(df_stat, peaks_info)


df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & abs(df_stat$Log2FC) > cutoff_Log2FC), ]
DARs_significant[[contrast_name]] <- df_stat_significant


df <- data.frame(seqnames = DARs_significant[[contrast_name]]$seqnames,
                 starts = DARs_significant[[contrast_name]]$start - 1,
                 ends = DARs_significant[[contrast_name]]$end, 
                 strands = DARs_significant[[contrast_name]]$strand)

write.table(df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)



############################ 2. Stem: POU1vsPP ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'Stem_POU1vsPP'
group.1 <- "Stem-POU1"
group.2 <- "Stem-PP"
annotation_group <- "anno1.loc"
cutoff_FDR <- 0.1 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 0.5
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

# export DARs in this contrast
df_stat <- archr_helper_markerPeaks_converter(DARs[[contrast_name]])
archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                   start(proj@peakSet),'-',end(proj@peakSet)))
idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
df_stat <- bind_cols(df_stat, peaks_info)


df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & abs(df_stat$Log2FC) > cutoff_Log2FC), ]
DARs_significant[[contrast_name]] <- df_stat_significant


df <- data.frame(seqnames = DARs_significant[[contrast_name]]$seqnames,
                 starts = DARs_significant[[contrast_name]]$start - 1,
                 ends = DARs_significant[[contrast_name]]$end, 
                 strands = DARs_significant[[contrast_name]]$strand)

write.table(df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)



############################ 3. Stem: POU2vsPP ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'Stem_POU2vsPP'
group.1 <- "Stem-POU2"
group.2 <- "Stem-PP"
annotation_group <- "anno1.loc"
cutoff_FDR <- 0.1 # cutoff is for peak selection, cutoffs for motifs are specificied in the text
cutoff_Log2FC <- 0.5
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)


# export DARs in this contrast
df_stat <- archr_helper_markerPeaks_converter(DARs[[contrast_name]])
archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
                                   start(proj@peakSet),'-',end(proj@peakSet)))
idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
df_stat <- bind_cols(df_stat, peaks_info)


df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & abs(df_stat$Log2FC) > cutoff_Log2FC), ]
DARs_significant[[contrast_name]] <- df_stat_significant


df <- data.frame(seqnames = DARs_significant[[contrast_name]]$seqnames,
                 starts = DARs_significant[[contrast_name]]$start - 1,
                 ends = DARs_significant[[contrast_name]]$end, 
                 strands = DARs_significant[[contrast_name]]$strand)

write.table(df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)












