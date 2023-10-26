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
output.dir <- '~/yuzhao1/work/final_RC2atac/peaks/4homer/output/'

# # read files
# proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
# DARs <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_EC.rds')
DARs_plot <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/EC_DAR_regions_DifferentContrastLists_ScriptInManuFolder.rds')


############################ star 1: colon core ###############################
contrast_name <- 'colon_core'
dir.create(paste0(output.dir, contrast_name))
# colon core is defined by the intersection of three groups: ACvsTI, POU2 vs PP, POU2vsPOU1

usefulpeaks_vector <- DARs_plot$EC_ACvsTI %>% 
  intersect(., DARs_plot$EC_POU2vsPP)%>%
  intersect(., DARs_plot$EC_POU2vsPOU1)

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)




############################ star 2: ileum core ###############################
contrast_name <- 'ileum_core'
dir.create(paste0(output.dir, contrast_name))
# ileum core is defined by the intersection of three groups: TIvsAC, PPvsPOU2, POU1vsPOU2

usefulpeaks_vector <- DARs_plot$EC_TIvsAC %>% 
  intersect(., DARs_plot$EC_PPvsPOU2)%>%
  intersect(., DARs_plot$EC_POU1vsPOU2)

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)



############################ star 3: colon only (intrinsic in colon, not gained in pouch2) ###############################
contrast_name <- 'colon_only'
dir.create(paste0(output.dir, contrast_name))
# defined by: ACvsTI minus POU2vsPP

usefulpeaks_vector <- DARs_plot$EC_ACvsTI %>% 
  setdiff(., DARs_plot$EC_POU2vsPP) # 45564 elements, too many

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)


############################ star 4: pouch2 only (gained in pouch from pp, but not colon features) ##########
contrast_name <- 'pouch2_only'
dir.create(paste0(output.dir, contrast_name))
# defined by: POU2vsTI minus ACvsTI

usefulpeaks_vector <- DARs_plot$EC_POU2vsPP %>% 
  setdiff(., DARs_plot$EC_ACvsTI)

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)

############################ star 5: ileum only (intrinsic in ileum, not preserved in pouch2) ###############################
contrast_name <- 'ileum_only'
dir.create(paste0(output.dir, contrast_name))
# defined by: TIvsAC minus PPvsPOU2

usefulpeaks_vector <- DARs_plot$EC_TIvsAC %>% 
  setdiff(., DARs_plot$EC_PPvsPOU2) # 26864 elements, too many

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)

########################### extra1: compare pouch2vsAC ###############################
contrast_name <- 'pouch2vsAC'
dir.create(paste0(output.dir, contrast_name))
usefulpeaks_vector <- DARs_plot$EC_POU2vsAC 

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)


########################### extra2: pouch2vsAC - TIvsAC ###############################
contrast_name <- 'pouch2vsAC_excludingTIvsAC'
dir.create(paste0(output.dir, contrast_name))

usefulpeaks_vector <- DARs_plot$EC_POU2vsAC %>% 
  setdiff(., DARs_plot$EC_TIvsAC) 

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)

########################### extra3: compare ACvspouch2 ###############################
contrast_name <- 'ACvspouch2'
dir.create(paste0(output.dir, contrast_name))
usefulpeaks_vector <- DARs_plot$EC_ACvsPOU2

usefulpeaks_df <- data.frame(seqnames = strsplit(usefulpeaks_vector, ':') %>% sapply(.,`[[`,1),
                             starts = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric() - 1,
                             ends = usefulpeaks_vector %>% strsplit(., ':') %>% sapply(.,`[[`,2) %>% strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                             strands = '*')
write.table(usefulpeaks_df, file=paste0(input.dir, contrast_name, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)






############################ below is time-consuming procedure ###############################

############################ 1. EC: ACvsTI ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'EC_ACvsTI'
group.1 <- "EC-AC"
group.2 <- "EC-TI"
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


############################ 2. EC: POU1vsPP ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'EC_POU1vsPP'
group.1 <- "EC-POU1"
group.2 <- "EC-PP"
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



############################ 3. EC: POU2vsPP ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'EC_POU2vsPP'
group.1 <- "EC-POU2"
group.2 <- "EC-PP"
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



############################ 4. EC: POU2vsPOU1 ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'EC_POU2vsPOU1'
group.1 <- "EC-POU2"
group.2 <- "EC-POU1"
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



############################ 1. EC: ACvsTI ###############################
temp_lineage <- 'epithelial'
contrast_name <- 'EC_ACvsTI'
group.1 <- "EC-AC"
group.2 <- "EC-TI"
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







