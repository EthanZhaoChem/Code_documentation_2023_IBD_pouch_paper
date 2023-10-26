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
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC/")

# proj <- addCoAccessibility(
#   ArchRProj = proj,
#   reducedDims = "Harmony",
#   dimsToUse = 1:99
# )

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "Harmony",
  knnIteration = 2000,
  dimsToUse = 1:99
)
  
saveArchRProject(ArchRProj = proj, load = T)










