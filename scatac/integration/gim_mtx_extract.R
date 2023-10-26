# filter again based on unconstrained predicted labels: at least they don't work based on gene score
dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(1)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'

# read files
proj.all <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj.epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj.tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj.bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj.myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj.others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")


gim <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneIntegrationMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

gim_imputed <- imputeMatrix(
  mat = gsm@assays@data$GeneIntegrationMatrix,
  imputeWeights = getImputeWeights(proj),
  threads = getArchRThreads(),
  verbose = FALSE,
  logFile = createLogFile("imputeMatrix")
)

saveRDS(gim, paste0(proj@projectMetadata$outputDirectory, '/GeneIntegrationMatrix.rds'))




