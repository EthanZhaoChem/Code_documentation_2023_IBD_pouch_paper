dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(limma)
library(edgeR)
library(variancePartition)
library(BiocParallel)
library(tidyverse)

source('~/yuzhao1/scripts/helper_archr.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')
addArchRThreads(4)
df_customized_enrichment <- list()

# 1. read files
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
DARs_plot <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/EC_DAR_regions_DifferentContrastLists_ScriptInManuFolder.rds')


# 2. analyze dar results: pp vs pou2
de_results <- readRDS('~/yuzhao1/work/final_RC2atac/dar_great/rds/da_results_POU_EC2vsPP.rds')

fclist <- topTable(
  de_results,
  coef = "anno1.subPP_EC",
  P.Value = 0.001,
  sort.by = 'logFC',
  number = 300000)

peaks_PPvsPOU2 <- rownames(fclist[which(fclist$adj.P.Val<0.05 & fclist$logFC > 1), ])
peaks_POU2vsPP <- rownames(fclist[which(fclist$adj.P.Val<0.05 & fclist$logFC < -1), ])

#
xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = peaks_PPvsPOU2)
df_customized_enrichment[['PPvsPOU2']] <- xx

#
xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = peaks_POU2vsPP)
df_customized_enrichment[['POU2vsPP']] <- xx


# 3. analyze dar results ac ti
de_results <- readRDS('~/yuzhao1/work/final_RC2atac/dar_great/rds/da_results_AC_ECvsTI_EC.rds')

fclist <- topTable(
  de_results,
  coef = "anno1.subTI_EC",
  sort.by = 'logFC',
  number = 300000)

peaks_TIvsAC <- rownames(fclist[which(fclist$adj.P.Val<0.05 & fclist$logFC > 1), ])
peaks_ACvsTI <- rownames(fclist[which(fclist$adj.P.Val<0.05 & fclist$logFC < -1), ])

#
xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = peaks_TIvsAC)
df_customized_enrichment[['TIvsAC']] <- xx

#
xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = peaks_ACvsTI)
df_customized_enrichment[['ACvsTI']] <- xx


# 4. analyze dar results pou1 pou2
de_results <- readRDS('~/yuzhao1/work/final_RC2atac/dar_great/rds/da_results_POU_EC2vsPOU_EC1.rds')

fclist <- topTable(
  de_results,
  coef = "anno1.subPOU_EC2",
  sort.by = 'logFC',
  number = 300000)

peaks_POU2vsPOU1 <- rownames(fclist[which(fclist$adj.P.Val<0.05 & fclist$logFC > 1), ])
peaks_POU1vsPOU2 <- rownames(fclist[which(fclist$adj.P.Val<0.05 & fclist$logFC < -1), ])

#
xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = peaks_POU2vsPOU1)
df_customized_enrichment[['POU2vsPOU1']] <- xx

#
xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = peaks_POU1vsPOU2)
df_customized_enrichment[['POU1vsPOU2']] <- xx


# 5. compare with wilcoxon
intersect(DARs_plot$EC_ACvsTI %>% gsub(':','_',.) %>% gsub('-','_',.), peaks_ACvsTI)
intersect(DARs_plot$EC_TIvsAC %>% gsub(':','_',.) %>% gsub('-','_',.), peaks_TIvsAC)


intersect(DARs_plot$EC_POU2vsPOU1 %>% gsub(':','_',.) %>% gsub('-','_',.), peaks_POU2vsPOU1)
intersect(DARs_plot$EC_POU1vsPOU2 %>% gsub(':','_',.) %>% gsub('-','_',.), peaks_POU1vsPOU2)











