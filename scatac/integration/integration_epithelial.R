# filter again based on unconstrained predicted labels: at least they don't work based on gene score

dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(12)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'

# read files
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2")
rna_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/epithelial3.rds')

## remove old labels
meta <- proj@cellColData
idx <- grep('predicted', colnames(meta))
proj@cellColData[, idx] <- NULL

# lineage constrained
df_rna <- rna_ref@meta.data
df_rna$embedding1 <- data.frame(rna_ref@reductions$harmony_umap@cell.embeddings)$UMAP_1
df_rna$embedding2 <- data.frame(rna_ref@reductions$harmony_umap@cell.embeddings)$UMAP_2
df_rna <- df_rna[, c('embedding1', 'embedding2')]

df_atac <- data.frame((proj@cellColData))
df_atac$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df_atac$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df_atac <- df_atac[, c('embedding1', 'embedding2')]


proj<- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = rna_ref,
  addToArrow = T,
  groupList = NULL,
  groupRNA = "anno1",
  nameCell = "predictedCell_lineage_constrained",
  nameGroup = "predictedGroup_lineage_constrained",
  nameScore = "predictedScore_lineage_constrained",
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  dimsToUse = 1:50,
  embeddingATAC = df_atac,
  embeddingRNA = df_rna,
  force = T
)


saveArchRProject(ArchRProj = proj, load = T)


# subset to EC after adding integrated gene score to arrow files

idxPass <- which(proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2"))
cellsPass <- proj$cellNames[idxPass]
proj.ec <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)


proj.ec <- addImputeWeights(proj_ec)
saveArchRProject(proj.ec, load = T)






















