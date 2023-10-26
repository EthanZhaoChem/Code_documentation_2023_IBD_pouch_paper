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
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem/")
proj <- addImputeWeights(proj)
saveArchRProject(ArchRProj = proj, load = T)


# ############################# Section0: save umap ################################
# save umap information to the proj folder

png('/project/gca/yuzhao1/work/final_RC2atac/annotation/plots/stem_ec/umap_anno2.png',width = 1800, height = 2000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno2
plot_df_umap_custom(df, show.label = 'name')
dev.off()

df <- data.frame(cell = proj$cellNames,
                 umap1 = proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`,
                 umap2 = proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`)
write.table(df, '~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem/harmony_umap.csv',sep = ',')



# ############################# Section1: Dim-Red ################################
# Dim Reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  dimsToUse = 2:100,
  varFeatures = 50000,
  clusterParams = list(
    sampleCells = 20000,
    resolution = c(2,2,2),
    n.start = 10
  ),
  sampleCellsPre = 20000,
  force = TRUE
)

saveArchRProject(ArchRProj = proj, load = T)



proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "patient",
  dimsToUse = 2:100,
  force = T
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "Harmony_UMAP",
  nNeighbors = 30,
  minDist = 0.3,
  force = T,
  dimsToUse = 1:30,
  metric = "cosine"
)

saveArchRProject(ArchRProj = proj, load = T)



# ############################# Section2: subset ################################
# in pouch
idxPass <- which(proj$biopsy_location %in% c('POU'))
cellsPass <- proj$cellNames[idxPass]
proj.EC_Stem_POU <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem_POU_version2",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T)
proj.EC_Stem_POU <- addImputeWeights(proj.EC_Stem_POU)
saveArchRProject(ArchRProj = proj.EC_Stem_POU, load = T)




# in AC TI
idxPass <- which(proj$biopsy_location %in% c('AC', 'TI'))
cellsPass <- proj$cellNames[idxPass]
proj.EC_Stem_ACTI <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem_ACTI_version2",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T)
proj.EC_Stem_ACTI <- addImputeWeights(proj.EC_Stem_ACTI)
saveArchRProject(ArchRProj = proj.EC_Stem_ACTI, load = T)


