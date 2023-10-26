# filter again based on unconstrained predicted labels: at least they don't work based on gene score

dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(4)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'


# read files

# proj.2 <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered")
# dic <- proj.2@cellColData
# dic <- data.frame(ID = rownames(dic),
#                   predictedGroup_Co_lineage=dic$predictedGroup_Co_lineage)
# proj$predictedGroup_Co_lineage <- mapvalues(proj$cellNames, from = dic$ID,
#                                             to = dic$predictedGroup_Co_lineage)

proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")



proj <- addImputeWeights(proj)
saveArchRProject(ArchRProj = proj, load = T)

# ############################# Section extra: anno1.loc update ################################
proj$anno1.loc <- paste0(proj$anno1, '-', proj$biopsy_location)

proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='POU'] <- 'Stromal-POU'
proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='PP'] <- 'Stromal-PP'
proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='TI'] <- 'Stromal-TI'
proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='AC'] <- 'Stromal-AC'

proj$anno1_loc <- proj$anno1.loc
saveArchRProject(ArchRProj = proj, load = T)


# # ############################# Section Last: anno1 ################################

df_annotation_res1.5 <- list(
  'C1' ='Glial',
  'C2' ='Pericyte',
  'C3' ='Myofibroblast',
  'C4' ='Stromal-1',
  'C5' ='Stromal-1',
  'C6' ='Stromal-2',
  'C7' ='Stromal-3',
  'C8' ='Lymphatic endothelium',
  'C9' ='Venous',
  'C10' ='Arterial',
  'C11' ='Arterial',
  'C12' ='Arterial',
  'C13' ='Venous',
  'C14' ='Arterial')

proj$anno1 <- mapvalues(proj$Harmony_Clusters_res1.5, from = names(df_annotation_res1.5), to = df_annotation_res1.5 %>%unlist())

png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_anno1.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_anno1_loc.png',width = 5000, height = 5200,res = 300)
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
    sampleCells = 3237,
    resolution = c(2,2),
    n.start = 10
  ),
  sampleCellsPre = 3237,
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
  minDist = 0.5,
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

png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_res1.5.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res1.5
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_res2.5.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res2.5
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_res5.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res5
plot_df_umap_custom(df, show.label = 'name')
dev.off()


# un-constrained labeling with epithelial rna cells only
rna_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/others3.rds')
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
  sampleCellsATAC = 3237,
  sampleCellsRNA = 4940,
  dimsToUse = 1:50,
  embeddingATAC = df_atac,
  embeddingRNA = df_rna,
  force = T
)
saveArchRProject(ArchRProj = proj, load = T)

png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_anno1_Co_lineage.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_Co_lineage
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_Co_lineage_loc.png',width = 5000, height = 5200,res = 300)
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
saveArchRProject(ArchRProj = proj, load = T)


png('~/yuzhao1/work/final_RC2atac/annotation/plots/others/umap_anno1_Co_Sample.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_Co_sample
plot_df_umap_custom(df, show.label = 'name')
dev.off()

# # ###################### Section 3: check markers #####################

gene <- c('IGHA2', 'IGHA1', 'CD27', 'SDC1', 'SLAMF7', 'CD3D', 'XBP1', 'PTPRC', 'LGR5', 'EPCAM',
          unlist(gca.epithelial.markers.anno2), unlist(gca.lineage.markers),
          unlist(gca.immune.markers.anno2), unlist(gca.mesenchymal.markers.anno2),
          gca.mesenchymal.markers %>% unlist())  %>% unique()

gene2 <- c('KLRF1', 'IFNG', 'FGFBP2', 'GZMB', 'CX3CR1', 'KLRG1', # NKT, has CD3D
           'IL7R', 'TCF7', 'PRKG1', 'PCDH9', 'AFF3', 'AREG','IL1R1', 'IL23R', 'KIT', # ILC  non-CD3D
           'GZMK', 'NCR1',  'KLRF1', # NK, non-CD3D, NON-IL7RN
           'SLC4A10', 'NCR3', 'KLRG1', 'GZMK', 'KLRG1') %>% unique() # MAIT, has CD3
gene <- c(gene, gene2)%>% intersect(., proj@geneAnnotation$genes$symbol)

gene_glial<- c('GRIK3',
                'MPZ',
                'NRXN1',
                'MAL',
                'AC016766.1',
                'SCN7A',
                'CADM2',
                'COL8A1',
                'SPP1',
                'TENM3',
                'MYOT',
                'GFRA3',
                'PPP2R2B',
                'RASGEF1C', 'S100B', 'NRXN1')%>% intersect(., proj@geneAnnotation$genes$symbol)

gene_others <- c('PECAM1', 'HEY1', 'EFNB2', 'ACKR1', 'VWF', 'ADAMDEC1', 'NOTCH3', 'MCAM',  'RGS5',
                 'ACTA2', 'TAGLN', 'DCN', 'DES', 'PROX1', 'LYVE1', 'CCL21', 'PLN', 'RERGL', 'KCNAB1', 'TAGLN', 'CNN1') 

gene <- unique(c(gene_glial, gene, gene_others))%>% intersect(., proj@geneAnnotation$genes$symbol)

p1 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = c(gene),
  embedding = "Harmony_UMAP",
  quantCut = c(0.01, 0.99),
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p1,
        name = paste0(c('Feature_umap', ".pdf")),
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
        name = paste0(c('Violin_res1.5', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

# # # # # # # ###################### Section 4: Remove low quality clusters #####################
# # filtered based on predicted label center (epithelial)
# # may need to take care of plasma cells later
# idxPass <- which(!proj$Harmony_Clusters_res1.5 %in% c('C2', 'C3'))
# cellsPass <- proj$cellNames[idxPass]
# proj.filtered <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = cellsPass,
#   outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2",
#   dropCells = TRUE,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = T
# )
# 


























