dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/helper_archr.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')
addArchRThreads(12)
# our.dir <- '~/yuzhao1/work/manu/rc2/plots/4atac_MarkerPeaks/'

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj_tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj_bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj_myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj_others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")
proj_immune <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_immune_filtered/")


for (temp_lineage in c("tcell", "bcell", "myeloid", "others", "immune")) {
  proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
  proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "Harmony",
    dimsToUse = 1:99
  )
  
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "Harmony",
    dimsToUse = 1:99
  )
  
  saveArchRProject(ArchRProj = proj, load = T)
}









