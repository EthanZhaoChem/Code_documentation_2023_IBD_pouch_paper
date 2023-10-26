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
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj <- addImputeWeights(proj)
saveArchRProject(ArchRProj = proj, load = T)

# ############################# Section extra: anno1.loc update ################################
proj$anno1.loc <- paste0(proj$anno1, '-', proj$biopsy_location)
proj$anno1_loc <- proj$anno1.loc
saveArchRProject(ArchRProj = proj, load = T)


# ############################# Section Last: anno1 ################################
df_annotation_res2 <- list(
  'C1' ='DC',
  'C2' ='Monocyte',
  'C3' ='Monocyte',
  'C4' ='Macrophage',
  'C5' ='Macrophage',
  'C6' ='Macrophage',
  'C7' ='Mast',
  'C8' ='Mast',
  'C9' ='Mast',
  'C10' ='Neutrophil',
  'C11' ='Neutrophil',
  'C12' ='Neutrophil')



proj$anno1 <- mapvalues(proj$Harmony_Clusters_res2, from = names(df_annotation_res2), to = df_annotation_res2 %>%unlist())

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_anno1.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_anno1_loc.png',width = 5000, height = 5200,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'name') +
  facet_wrap(~ biopsy_location) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
dev.off()
saveArchRProject(ArchRProj = proj, load = T)

############################# Section1: Dim-Red ################################
# Dim Reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  dimsToUse = 2:50,
  varFeatures = 50000,
  clusterParams = list(
    sampleCells = 3000,
    resolution = c(2, 2),
    n.start = 10
  ),
  sampleCellsPre = 3000,
  force = TRUE
)

saveArchRProject(ArchRProj = proj, load = T)



proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "patient",
  dimsToUse = 2:50,
  force = T
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "Harmony_UMAP",
  nNeighbors = 50,
  minDist = 0.05,
  force = T,
  dimsToUse = 1:49,
  metric = "cosine"
)



proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_res1.5",
  force = T,
  dimsToUse = 1:49,
  resolution = 1.5,
  maxClusters = 50,
)

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_res2",
  force = T,
  dimsToUse = 1:49,
  resolution = 2,
  maxClusters = 50,
)

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_res2.5",
  force = T,
  dimsToUse = 1:49,
  resolution = 2.5,
  maxClusters = 50,
)
proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_res5",
  force = T,
  dimsToUse = 1:49,
  resolution = 5,
  maxClusters = 50,
)
saveArchRProject(ArchRProj = proj, load = T)

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_res1.5.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res1.5
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_res2.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res2
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_res2.5.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res2.5
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_res5.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res5
plot_df_umap_custom(df, show.label = 'name')
dev.off()


# un-constrained labeling with epithelial rna cells only
rna_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/myeloid3.rds')
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
  nameCell = "predictedCell_Co_lineage",
  nameGroup = "predictedGroup_Co_lineage",
  nameScore = "predictedScore_Co_lineage",
  sampleCellsATAC = 3000,
  sampleCellsRNA = 3000,
  dimsToUse = 1:50,
  embeddingATAC = df_atac,
  embeddingRNA = df_rna,
  force = T
)



saveArchRProject(ArchRProj = proj, load = T)

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_anno1_Co_lineage.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_Co_lineage
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/myeloid/umap_anno1_Co_lineage_loc.png',width = 5000, height = 5200,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_Co_lineage
plot_df_umap_custom(df, show.label = 'name') +
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

gene2 <- c('KLRF1', 'IFNG', 'FGFBP2', 'GZMB', 'CX3CR1', 'KLRG1', # NKT, has CD3D
           'IL7R', 'TCF7', 'PRKG1', 'PCDH9', 'AFF3', 'AREG','IL1R1', 'IL23R', 'KIT', # ILC  non-CD3D
           'GZMK', 'NCR1',  'KLRF1', # NK, non-CD3D, NON-IL7RN
           'SLC4A10', 'NCR3', 'KLRG1', 'GZMK', 'KLRG1') %>% unique() # MAIT, has CD3
gene <- c(gene, gene2)%>% intersect(., proj@geneAnnotation$genes$symbol)

gene_cDC1 <- c('AUTS2', 'HDAC9', 'IDO1', 'CCSER1', 'WDFY4','CADM1',
          'RGCC', 'CAMK2D', 'NEGR1', 'CPNE3','CLEC9A')%>% intersect(., proj@geneAnnotation$genes$symbol)

gene_neutrophil <- c('G0S2', 'CXCL8', 'CSF3R', 'S100A9', 'S100A8',
                    'IFITM2', 'FCGR3B', 'LUCAT1', 'AQP9', 'CXCR4',
                    'SLC25A37', 'IL1R2', 'LIMK2', 'NIBAN1', 'EPHB1')%>% intersect(., proj@geneAnnotation$genes$symbol)

gene_mast <- c('TPSAB1', 'TPSB2', 'CPA3', 'CTSG', 'HPGDS', 'SLC24A3',
               'KIT', 'CD69', 'MS4A2', 'HPGD', 'CADPS', 'MAST4', 'LTC4S', 'BACE2')

gene_LymDC <- c('CXCL9', 'TXN','CRIP1', 'MARCKSL1', 'CCL19','IDO1','LAMP3',
                'FSCN1','CCR7','PPA1','EBI3')



p1 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = c(gene),
  embedding = "Harmony_UMAP",
  quantCut = c(0.01, 0.99),
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p1,
        name = paste0(c('Feature_umap2_', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Harmony_Clusters_res1.5",
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
        name = paste0(c('Violin_res1.5_2_', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

# # # # # # # ###################### Section 4: Remove low quality clusters #####################
# # filtered based on predicted label center (epithelial)
# # may need to take care of plasma cells later
# idxPass <- which(!proj$Harmony_Clusters_res1.5 %in% c('C10'))
# cellsPass <- proj$cellNames[idxPass]
# proj.filtered <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = cellsPass,
#   outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2",
#   dropCells = TRUE,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = T
# )



























