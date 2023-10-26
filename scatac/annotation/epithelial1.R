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
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial/")
# proj <- addImputeWeights(proj)
# saveArchRProject(ArchRProj = proj, load = T)
# ############################# Section1: Dim-Red ################################
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
# proj <- addUMAP(
#   ArchRProj = proj,
#   reducedDims = "Harmony",
#   name = "Harmony_UMAP",
#   nNeighbors = 30,
#   minDist = 0.2,
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
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_Co_sample.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$predictedGroup_Co_sample
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_loc_Co_sample.png',width = 5000, height = 5200,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$predictedGroup_Co_sample
# plot_df_umap_custom(df, show.label = 'name')+
#     facet_wrap(~ biopsy_location) +
#     theme(
#       strip.background = element_rect(fill = "white", colour = "white"),
#       strip.text = element_text(size = 12)
#     )
# dev.off()
#
# # un-constrained labeling
# rna_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
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
#   nameCell = "predictedCell_un",
#   nameGroup = "predictedGroup_un",
#   nameScore = "predictedScore_un",
#   sampleCellsATAC = 10000,
#   sampleCellsRNA = 10000,
#   dimsToUse = 1:50,
#   embeddingATAC = df_atac,
#   embeddingRNA = df_rna,
#   force = T
# )
# saveArchRProject(ArchRProj = proj, load = T)
# 
# png('~/yuzhao1/work/final_RC2atac/annotation/plots/epithelial/umap_un.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$predictedGroup_un
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()

################# Section 2: Sample constrained Integration ###########
# rna_ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
# rna_ref$Sample_ID_Corrected <- paste0(rna_ref$Patient_ID,'-', rna_ref$biopsy_location)
# rna_cellnames_per_sample <- list()
# atac_cellnames_per_sample <- list()
# pairs <- list()
# 
# # check whether sample names are shared between atac and rna
# sample_names <- unique(rna_ref$Sample_ID_Corrected)
# sum(unique(rna_ref$Sample_ID_Corrected) %in% unique(proj$Sample_corrected))
# sum(unique(proj$Sample_corrected) %in% unique(rna_ref$Sample_ID_Corrected))
# 
# 
# for (sample.temp in sample_names){
#   id.rna <- which(rna_ref$Sample_ID_Corrected == sample.temp)
#   rna_cellnames_per_sample[[sample.temp]] <- Cells(rna_ref)[id.rna]
#   
#   id.atac <- which(proj$Sample_corrected == sample.temp)
#   atac_cellnames_per_sample[[sample.temp]] <- proj$cellNames[id.atac]
#   
#   pairs[[sample.temp]] <- SimpleList(ATAC = atac_cellnames_per_sample[[sample.temp]],
#                                      RNA = rna_cellnames_per_sample[[sample.temp]])
# }
# 
# # add list
# groupList <- SimpleList(pairs)
# 
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
#   groupList = groupList,
#   groupRNA = "anno1",
#   nameCell = "predictedCell_Co_sample",
#   nameGroup = "predictedGroup_Co_sample",
#   nameScore = "predictedScore_Co_sample",
#   sampleCellsATAC = 10000,
#   sampleCellsRNA = 10000,
#   dimsToUse = 1:50,
#   embeddingATAC = df_atac,
#   embeddingRNA = df_rna,
#   force = T
# )

# saveArchRProject(ArchRProj = proj, load = T)
# 
# png('~/yuzhao1/work/final_RC2atac/preprocess/plots/umap_Co_anno1_loc.png',width = 5000, height = 5200,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$predictedGroup_Co_sample
# plot_df_umap_custom(df, show.label = 'name')+
#   facet_wrap(~ biopsy_location) +
#   theme(
#     strip.background = element_rect(fill = "white", colour = "white"),
#     strip.text = element_text(size = 12)
#   )
# dev.off()
# 
# 
# png('~/yuzhao1/work/final_RC2atac/preprocess/plots/umap_Co_lineage.png',width = 2600, height = 3000,res = 300)
# df <- data.frame(proj@cellColData)
# df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
# df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# df$cluster_name <- proj$predictedGroup_lineage_Co_sample
# plot_df_umap_custom(df, show.label = 'name')
# dev.off()


###################### Section3: Check clusters #####################
# all.gsm <- getMatrixFromProject(
#   ArchRProj = proj,
#   useMatrix = "GeneScoreMatrix",
#   useSeqnames = NULL,
#   verbose = TRUE,
#   binarize = FALSE,
#   threads = 1,
#   logFile = createLogFile("getMatrixFromProject")
# )
# 
# all.gsm.imputed <- imputeMatrix(
#   mat = all.gsm@assays@data$GeneScoreMatrix,
#   imputeWeights = getImputeWeights(proj),
#   threads = 1,
#   verbose = T,
#   logFile = createLogFile("imputeMatrix")
# )
# 
# library(Matrix)
# dgc_imput_mat <- as(all.gsm.imputed, "dgCMatrix")
# # all.gsm@assays@data$ImputedGSM <- dgc_imput_mat
# # only save imputed gene score matrix
# saveRDS(dgc_imput_mat, paste0(proj@projectMetadata$outputDirectory, '/gsm_dgc_imput_mat.rds'))


# markersGS <- getMarkerFeatures(
#   ArchRProj = proj, 
#   useMatrix = "GeneScoreMatrix", 
#   groupBy = "Harmony_Clusters_res1.5",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon",
#   useGroups = c('C19')
#   #,bgdGroups = c('')
# )
# 
# markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
# markerList$C19
# 
# 
# 
# 


gene <- c('IGHA2', 'IGHA1', 'CD27', 'SDC1', 'SLAMF7', 'CD3D', 'XBP1', 'PTPRC', 'LGR5', 'EPCAM', unlist(gca.epithelial.markers), unlist(gca.lineage.markers))  %>% unique()
gene <- intersect(gene, proj@geneAnnotation$genes$symbol)

gene <- gca.mesenchymal.markers %>% unlist()

p1 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = c(gene),
  embedding = "Harmony_UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(proj)
  )

plotPDF(plotList = p1,
        name = paste0(c('Feature_umap2', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Harmony_Clusters_res2.5",
  colorBy = "GeneScoreMatrix",
  name = gene,
  plotAs = "violin",
  alpha = 0.4,
  maxCells = 10000,
  addBoxPlot = TRUE,
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p2,
        name = paste0(c('Violin2', ".pdf")),
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)



###################### Section 4: Remove low quality clusters #####################
# filtered based on predicted label center (immune and stromal)
# may need to take care of plasma cells later
idxPass <- which(!proj$Harmony_Clusters_res2.5 %in% c('C27', 'C31', 'C32'))
cellsPass <- proj$cellNames[idxPass]
proj.filtered <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass,
  outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
















