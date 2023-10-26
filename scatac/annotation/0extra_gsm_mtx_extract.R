# filter again based on unconstrained predicted labels: at least they don't work based on gene score
dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
library(Matrix)
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

for (temp_lineage in c('all', 'epithelial', 'tcell', 'bcell', 'myeloid', 'others')) {
  proj <- paste0('proj.', temp_lineage) %>% as.name(.) %>% eval(.)
  cat(paste0('Geting GSM: ', temp_lineage, '\n'))
  gsm <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
  )
  
  cat(paste0('Imputing GSM: ', temp_lineage, '\n'))
  gsm_imputed <- imputeMatrix(
    mat = gsm@assays@data$GeneScoreMatrix,
    imputeWeights = getImputeWeights(proj),
    threads = getArchRThreads(),
    verbose = T,
    logFile = createLogFile("imputeMatrix")
  )
  
  cat(paste0('Saving imputed GSM_imputed: ', temp_lineage, '\n'))
  rownames(gsm_imputed) <- gsm@elementMetadata$name
  saveRDS(gsm_imputed, paste0(proj@projectMetadata$outputDirectory, '/GeneScoreMatrix_imputed.rds'))
  
  gc()
}





