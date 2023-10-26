dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/helper_archr.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")

addArchRThreads(4)
out.dir <- '/project/gca/yuzhao1/work/final_RC2atac/peaks/customized_enrichment/'

# # read files
# DARs <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/3DARs/DARs_EC.rds')
DARs_plot <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/EC_DAR_regions_DifferentContrastLists_ScriptInManuFolder.rds')
df_customized_enrichment <- list()

############################ star 1: colon core ###############################
contrast_name <- 'colon_core'
usefulpeaks_vector <- DARs_plot$EC_ACvsTI %>% 
  intersect(., DARs_plot$EC_POU2vsPP)%>%
  intersect(., DARs_plot$EC_POU2vsPOU1)

xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)

df_customized_enrichment[[contrast_name]] <- xx

############################ star 2: ileum core ###############################
contrast_name <- 'ileum_core'
usefulpeaks_vector <- DARs_plot$EC_TIvsAC %>% 
  intersect(., DARs_plot$EC_PPvsPOU2)%>%
  intersect(., DARs_plot$EC_POU1vsPOU2)


xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)

df_customized_enrichment[[contrast_name]] <- xx

############################ star 3: colon only (intrinsic in colon, not gained in pouch2) ###############################
contrast_name <- 'colon_only'
usefulpeaks_vector <- DARs_plot$EC_ACvsTI %>% 
  setdiff(., DARs_plot$EC_POU2vsPP) # 45564 elements, too many

xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)


df_customized_enrichment[[contrast_name]] <- xx

############################ star 4: pouch2 only (gained in pouch from pp, but not colon features) ##########
contrast_name <- 'pouch2_only'
usefulpeaks_vector <- DARs_plot$EC_POU2vsPP %>% 
  setdiff(., DARs_plot$EC_ACvsTI)

xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)

df_customized_enrichment[[contrast_name]] <- xx

############################ star 5: ileum only (intrinsic in ileum, not preserved in pouch2) ###############################
contrast_name <- 'ileum_only'
usefulpeaks_vector <- DARs_plot$EC_TIvsAC %>% 
  setdiff(., DARs_plot$EC_PPvsPOU2) # 26864 elements, too many

xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)
df_customized_enrichment[[contrast_name]] <- xx

########################### extra1: compare pouch2vsAC ###############################
contrast_name <- 'pouch2vsAC'
usefulpeaks_vector <- DARs_plot$EC_POU2vsAC 

xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)
df_customized_enrichment[[contrast_name]] <- xx


########################### extra2: pouch2vsAC - TIvsAC ###############################
contrast_name <- 'pouch2vsAC_excludingTIvsAC'
usefulpeaks_vector <- DARs_plot$EC_POU2vsAC %>% 
  setdiff(., DARs_plot$EC_TIvsAC) 

xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)

df_customized_enrichment[[contrast_name]] <- xx

########################### extra3: compare ACvspouch2 ###############################
contrast_name <- 'ACvspouch2'
usefulpeaks_vector <- DARs_plot$EC_ACvsPOU2

xx <- archr_customized_motif_enrichment(proj, 
                                        peakAnnotation = "Motif", 
                                        candidate_peaks_vector = usefulpeaks_vector)

df_customized_enrichment[[contrast_name]] <- xx


########################### extra all: plot all names in DAR ###############################
for (contrast_name in names(DARs_plot)){
  usefulpeaks_vector <- DARs_plot[[contrast_name]]
  xx <- archr_customized_motif_enrichment(proj, 
                                          peakAnnotation = "Motif", 
                                          candidate_peaks_vector = usefulpeaks_vector)
  
  df_customized_enrichment[[contrast_name]] <- xx
}


saveRDS(df_customized_enrichment, "/project/gca/yuzhao1/work/final_RC2atac/peaks/customized_enrichment/EC_allPossible_enrichments.rds")









