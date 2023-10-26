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
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered/")
# proj <- addImputeWeights(proj)
# saveArchRProject(ArchRProj = proj, load = T)
############################# Section1: Dim-Red ################################
# # Dim Reduction
# proj <- addIterativeLSI(
#   ArchRProj = proj,
#   useMatrix = "TileMatrix",
#   name = "IterativeLSI",
#   iterations = 6,
#   dimsToUse = 2:100,
#   varFeatures = 50000,
#   clusterParams = list(
#     sampleCells = 50000,
#     resolution = c(0.5, 1, 1.5, 2, 2),
#     n.start = 10
#   ),
#   sampleCellsPre = 50000,
#   force = TRUE
# )
# 
# saveArchRProject(ArchRProj = proj, load = T)
# 
# 
# 
# proj <- addHarmony(
#   ArchRProj = proj,
#   reducedDims = "IterativeLSI",
#   name = "Harmony",
#   groupBy = "patient",
#   dimsToUse = 2:100,
#   force = T
# )
# 
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "Harmony_UMAP",
  nNeighbors = 20,
  minDist = 0.2,
  force = T,
  dimsToUse = 1:99,
  metric = "cosine"
)

# 
# 
# proj <- addClusters(
#   input = proj,
#   reducedDims = "Harmony",
#   method = "Seurat",
#   name = "Harmony_Clusters_res1.5",
#   force = T,
#   dimsToUse = 1:99,
#   resolution = 1.5,
#   maxClusters = 100,
# )
# 
# proj <- addClusters(
#   input = proj,
#   reducedDims = "Harmony",
#   method = "Seurat",
#   name = "Harmony_Clusters_res2.5",
#   force = T,
#   dimsToUse = 1:99,
#   resolution = 2.5,
#   maxClusters = 100,
# )
# proj <- addClusters(
#   input = proj,
#   reducedDims = "Harmony",
#   method = "Seurat",
#   name = "Harmony_Clusters_res5",
#   force = T,
#   dimsToUse = 1:99,
#   resolution = 5,
#   maxClusters = 100,
# )
# saveArchRProject(ArchRProj = proj, load = T)
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_res1.5.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$Harmony_Clusters_res1.5
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_res2.5.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$Harmony_Clusters_res2.5
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_res5.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$Harmony_Clusters_res5
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()


# un-constrained labeling
rna_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/epithelial3.rds')
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
  addToArrow = F,
  groupList = NULL,
  groupRNA = "anno1",
  nameCell = "predictedCell_un",
  nameGroup = "predictedGroup_un",
  nameScore = "predictedScore_un",
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  dimsToUse = 1:50,
  embeddingATAC = df_atac,
  embeddingRNA = df_rna,
  force = T
)

proj<- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = rna_ref,
  addToArrow = F,
  groupList = NULL,
  groupRNA = "lineage",
  nameCell = "predictedCell_lineage_un",
  nameGroup = "predictedGroup_lineage_un",
  nameScore = "predictedScore_lineage_un",
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  dimsToUse = 1:50,
  embeddingATAC = df_atac,
  embeddingRNA = df_rna,
  force = T
)

saveArchRProject(ArchRProj = proj, load = T)

png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_un_loc.png',width = 5000, height = 5200,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_un
plot_df_umap_custom(df, show.label = 'name')+
  facet_wrap(~ biopsy_location) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
dev.off()

# ###################### Section 3: check markers #####################

gene <- c('IGHA2', 'IGHA1', 'CD27', 'SDC1', 'SLAMF7', 'CD3D', 'XBP1', 'PTPRC', 'LGR5', 'EPCAM',
          unlist(gca.epithelial.markers.anno2), unlist(gca.lineage.markers),
          unlist(gca.immune.markers.anno2), unlist(gca.mesenchymal.markers.anno2),
          gca.mesenchymal.markers %>% unlist())  %>% unique()
gene <- intersect(gene, proj@geneAnnotation$genes$symbol)


p1 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = c(gene),
  embedding = "Harmony_UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p1,
        name = paste0(c('Feature_umap', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "predictedGroup_un",
  colorBy = "GeneScoreMatrix",
  name = gene,
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 4,
  maxCells = 10000,
  addBoxPlot = TRUE,
  ratioYX = 0.6,
  imputeWeights = getImputeWeights(proj)
)


plotPDF(plotList = p2,
        name = paste0(c('Violin', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

pp <- plotGroups(
  ArchRProj = proj,
  groupBy = "predictedGroup_un",
  colorBy = "cellColData",
  name = c('nFrags','ReadsInTSS'),
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 4,
  maxCells = 10000,
  addBoxPlot = TRUE,
  ratioYX = 0.6,
  imputeWeights = getImputeWeights(proj)
)


plotPDF(plotList = pp,
        name = paste0(c('Violin_nFrags', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

# ###################### Section 4: Remove low quality clusters #####################
# filtered based on predicted label center (immune and stromal)
# may need to take care of plasma cells later
idxPass <- which(proj$predictedGroup_lineage_un %in% c('epithelial'))
cellsPass <- proj$cellNames[idxPass]
proj.filtered <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
















