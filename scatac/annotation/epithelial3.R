# filter again based on unconstrained predicted labels: at least they don't work based on gene score

dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(6)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'



# read files
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2")


# table(proj@cellColData[, c('Harmony_Clusters_res1.5', 'biopsy_location')])
# 
# # proj <- addImputeWeights(proj)
# # saveArchRProject(ArchRProj = proj, load = T)


# # ############################# Section extra: anno2 update ################################
# proj_stem <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_Stem/")
# cells1 <- proj_stem$cellNames[proj_stem$anno2 == 'Stem1'] 
# cells2 <- proj_stem$cellNames[proj_stem$anno2 == 'Stem2'] 
# 
# proj$anno2 <- proj$anno1
# idx1 <- match(cells1, proj$cellNames)
# proj$anno2[idx1] <- 'Stem1' 
# idx2 <- match(cells2, proj$cellNames)
# proj$anno2[idx2] <- 'Stem2' 
# 
# saveArchRProject(ArchRProj = proj, load = T)

 
# # ############################# Section extra: anno1.loc & anno1_loc update ################################
# proj$anno1.loc <- paste0(proj$anno1, '-', proj$biopsy_location)
# proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2") & proj$biopsy_location =='POU'] <- 'EC-POU1'
# proj$anno1.loc[proj$anno1 %in% c("EC2-1", "EC2-2") & proj$biopsy_location =='POU'] <- 'EC-POU2'
# proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2") & proj$biopsy_location =='PP'] <- 'EC-PP'
# proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2") & proj$biopsy_location =='TI'] <- 'EC-TI'
# proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2") & proj$biopsy_location =='AC'] <- 'EC-AC'
# 
# proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='POU'] <- 'Goblet-POU'
# proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='PP'] <- 'Goblet-PP'
# proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='TI'] <- 'Goblet-TI'
# proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='AC'] <- 'Goblet-AC'
# 
# proj$anno1.loc[proj$anno2 %in% c("Stem1") & proj$biopsy_location =='POU'] <- 'Stem-POU1'
# proj$anno1.loc[proj$anno2 %in% c("Stem2") & proj$biopsy_location =='POU'] <- 'Stem-POU2'
# 
# proj$anno1_loc <- proj$anno1.loc
# saveArchRProject(ArchRProj = proj, load = T)




# # ############################# Section Last: anno1 ################################
# df_annotation_res1.5 <- list(
#   'C1' ='Tuft',
#   'C2' ='EC1-2',
#   'C3' ='EC2-2',
#   'C4' ='EC1-1',
#   'C5' ='EC1-2',
#   'C6' ='EC1-2',
#   'C7' ='EC1-2',
#   'C8' ='EC1-2',
#   'C9' ='EC1-2',
#   'C10' ='EC1-2',
#   'C11' ='EEC',
#   'C12' ='Paneth',
#   'C13' ='Goblet1',
#   'C14' ='Goblet2',
#   'C15' ='Goblet2',
#   'C16' ='EC2-1',
#   'C17' ='BEST4',
#   'C18' ='EC2-1',
#   'C19' ='EC2-1',
#   'C20' ='EC2-1',
#   'C21' ='EC2-1',
#   'C22' ='EC2-1',
#   'C23' ='EC1-1',
#   'C24' ='EC1-1',
#   'C25' ='EC2-1',
#   'C26' ='Stem')
# proj$anno1 <- mapvalues(proj$Harmony_Clusters_res1.5, from = names(df_annotation_res1.5), to = df_annotation_res1.5 %>%unlist())
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_anno1.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$anno1
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_anno1_loc.png',width = 5000, height = 5200,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$anno1
# plot_df_umap_custom(df, show.label = 'name') +
#   facet_wrap(~ biopsy_location) +
#   theme(
#     strip.background = element_rect(fill = "white", colour = "white"),
#     strip.text = element_text(size = 12)
#   )
# dev.off()
# saveArchRProject(ArchRProj = proj, load = T)

# ############################# Section1: Dim-Red ################################
# # Dim Reduction
# proj <- addIterativeLSI(
#   ArchRProj = proj,
#   useMatrix = "TileMatrix",
#   name = "IterativeLSI",
#   iterations = 5,
#   dimsToUse = 2:100,
#   varFeatures = 150000,
#   clusterParams = list(
#     sampleCells = 150000,
#     resolution = c(2,2,2,2),
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
# proj <- addUMAP(
#   ArchRProj = proj,
#   reducedDims = "Harmony",
#   name = "Harmony_UMAP",
#   nNeighbors = 30,
#   minDist = 0.3,
#   force = T,
#   dimsToUse = 1:99,
#   metric = "cosine"
# )
# 
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
# 
# 
# # un-constrained labeling with epithelial rna cells only
# rna_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/epithelial3.rds')
# df_rna <- rna_ref@meta.data
# df_rna$embedding1 <- data.frame(rna_ref@reductions$harmony_umap@cell.embeddings)$UMAP_1
# df_rna$embedding2 <- data.frame(rna_ref@reductions$harmony_umap@cell.embeddings)$UMAP_2
# df_rna <- df_rna[, c('embedding1', 'embedding2')]
# 
# df_atac <- data.frame((proj@cellColData))
# df_atac$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df_atac$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df_atac <- df_atac[, c('embedding1', 'embedding2')]
# 
# proj<- addGeneIntegrationMatrix(
#   ArchRProj = proj,
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "Harmony",
#   seRNA = rna_ref,
#   addToArrow = F,
#   groupList = NULL,
#   groupRNA = "anno1",
#   nameCell = "predictedCell_Co_lineage",
#   nameGroup = "predictedGroup_Co_lineage",
#   nameScore = "predictedScore_Co_lineage",
#   sampleCellsATAC = 10000,
#   sampleCellsRNA = 10000,
#   dimsToUse = 1:50,
#   embeddingATAC = df_atac,
#   embeddingRNA = df_rna,
#   force = T
# )
# 
# 
# saveArchRProject(ArchRProj = proj, load = T)
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_Co_lineage.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$predictedGroup_Co_lineage
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_Co_lineage_loc.png',width = 5000, height = 5200,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$predictedGroup_Co_lineage
# plot_df_umap_custom(df, show.label = 'name') +
#   facet_wrap(~ biopsy_location) +
#   theme(
#     strip.background = element_rect(fill = "white", colour = "white"),
#     strip.text = element_text(size = 12)
#   )
# dev.off()
# 
# 
# # ###################### Section 3: check markers #####################
# 
# gene <- c('IGHA2', 'IGHA1', 'CD27', 'SDC1', 'SLAMF7', 'CD3D', 'XBP1', 'PTPRC', 'LGR5', 'EPCAM',
#           unlist(gca.epithelial.markers.anno2), unlist(gca.lineage.markers),
#           unlist(gca.immune.markers.anno2), unlist(gca.mesenchymal.markers.anno2),
#           gca.mesenchymal.markers %>% unlist())  %>% unique()
# gene <- intersect(gene, proj@geneAnnotation$genes$symbol)
# 
# gene_EC1.23 <- c('PLCG2', 'CD55', 'TMCC3', 'FHL2', 'KLF6',
#                  'REG1B',
#                  'REG1A',
#                  'CCL20',
#                  'CLCA4',
#                  'AC083837.1',
#                  'LCN2',
#                  'DUOX2',
#                  'NOS2',
#                  'BIRC3',
#                  'DUOXA2',
#                  'CEACAM6',
#                  'CASP5',
#                  'GCNT3',
#                  'NFKBIA',
#                  'NAALADL2',
#                  'RBFOX1',
#                  'PTPRK', 
#                  'FHIT')
# gene <- intersect(gene_EC1.23, proj@geneAnnotation$genes$symbol)
# 
# p1 <- plotEmbedding(
#   ArchRProj = proj,
#   colorBy = "GeneScoreMatrix",
#   name = c(gene),
#   embedding = "Harmony_UMAP",
#   quantCut = c(0.01, 0.95),
#   imputeWeights = getImputeWeights(proj)
# )
# 
# plotPDF(plotList = p1,
#         name = paste0(c('Feature_umap_2', ".pdf")),
#         ArchRProj = proj,
#         addDOC = FALSE, width = 5, height = 5)
# 
# p2 <- plotGroups(
#   ArchRProj = proj,
#   groupBy = "Harmony_Clusters_res1.5",
#   colorBy = "GeneScoreMatrix",
#   name = gene,
#   plotAs = "violin",
#   alpha = 0.4,
#   baseSize = 4,
#   maxCells = 10000,
#   addBoxPlot = TRUE,
#   ratioYX = 0.6,
#   imputeWeights = getImputeWeights(proj)
# )
# 
# 
# plotPDF(plotList = p2,
#         name = paste0(c('Violin_res1.5_2', ".pdf")),
#         ArchRProj = proj,
#         addDOC = FALSE, width = 5, height = 5)

# # ###################### Section 4: Remove low quality clusters #####################
# # filtered based on predicted label center (immune and stromal)
# # may need to take care of plasma cells later
# idxPass <- which(proj$predictedGroup_lineage_un %in% c('epithelial'))
# cellsPass <- proj$cellNames[idxPass]
# proj.filtered <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = cellsPass,
#   outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2",
#   dropCells = TRUE,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = T
# )

# # ###################### Section 5: increase stem resolution #####################

# idxPass <- which(proj$anno1 %in% c('Stem'))
# cellsPass <- proj$cellNames[idxPass]
# proj.stem <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = cellsPass,
#   outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_Stem",
#   dropCells = TRUE,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = T
# )



idxPass <- which(proj$anno1 %in% c("Goblet1", "Goblet2"))
cellsPass <- proj$cellNames[idxPass]
proj.Goblet <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_Goblet",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.Goblet <- addImputeWeights(proj.Goblet)
saveArchRProject(ArchRProj = proj.Goblet, load = T)




idxPass <- which(proj$anno1 %in% c("Stem", "EC1-1", "EC1-2", "EC2-1", "EC2-2"))
cellsPass <- proj$cellNames[idxPass]
proj.EC_Stem <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.EC_Stem <- addImputeWeights(proj.EC_Stem)
saveArchRProject(ArchRProj = proj.EC_Stem, load = T)


# EC and Stem in AC TI
idxPass <- which(proj$anno1 %in% c("Stem", "EC1-1", "EC1-2", "EC2-1", "EC2-2")   &
                   proj$biopsy_location %in% c('AC', 'TI'))
cellsPass <- proj$cellNames[idxPass]
proj.EC_Stem_ACTI <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem_ACTI",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.EC_Stem_ACTI <- addImputeWeights(proj.EC_Stem_ACTI)
saveArchRProject(ArchRProj = proj.EC_Stem_ACTI, load = T)


# EC and Stem in POU PP
idxPass <- which(proj$anno1 %in% c("Stem", "EC1-1", "EC1-2", "EC2-1", "EC2-2")   &
                   proj$biopsy_location %in% c('PP', 'POU'))
cellsPass <- proj$cellNames[idxPass]
proj.EC_Stem_PPPOU <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem_PPPOU",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.EC_Stem_PPPOU <- addImputeWeights(proj.EC_Stem_PPPOU)
saveArchRProject(ArchRProj = proj.EC_Stem_PPPOU, load = T)



# EC and Stem in POU
idxPass <- which(proj$anno1 %in% c("Stem", "EC1-1", "EC1-2", "EC2-1", "EC2-2")   &
                   proj$biopsy_location %in% c('POU'))
cellsPass <- proj$cellNames[idxPass]
proj.EC_Stem_POU <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC_Stem_POU",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.EC_Stem_POU <- addImputeWeights(proj.EC_Stem_POU)
saveArchRProject(ArchRProj = proj.EC_Stem_POU, load = T)







