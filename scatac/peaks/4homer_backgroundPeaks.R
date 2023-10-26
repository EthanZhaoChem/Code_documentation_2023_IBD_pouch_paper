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

# create background peak set for each lineage
for (temp_lineage in  unique(proj_union$lineage)) {
  proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
  peaks.all <- proj@peakSet
  df <- data.frame(seqnames = seqnames(peaks.all),
                   starts = start(peaks.all) - 1,
                   ends = end(peaks.all), 
                   strands = strand(peaks.all))
  
  write.table(df, file=paste0('~/yuzhao1/work/final_RC2atac/peaks/4homer/backgroundPeaks/',
                              temp_lineage, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
}











