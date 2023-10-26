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
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_Stem/")

# table(proj@cellColData[, c('Harmony_Clusters_res1.5', 'biopsy_location')])
# 
# proj <- addImputeWeights(proj)
# saveArchRProject(ArchRProj = proj, load = T)


# ############################# Section Last: anno1 ################################
# df_annotation_res1.5 <- list()
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
# Dim Reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  dimsToUse = 2:100,
  varFeatures = 50000,
  clusterParams = list(
    sampleCells = 1702,
    resolution = c(2,2),
    n.start = 10
  ),
  sampleCellsPre = 1702,
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

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_res1.5",
  force = T,
  dimsToUse = 1:30,
  resolution = 1.5,
  maxClusters = 100,
)

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_res0.5",
  force = T,
  dimsToUse = 1:30,
  resolution = 0.5,
  maxClusters = 100,
)

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_res1",
  force = T,
  dimsToUse = 1:30,
  resolution = 1,
  maxClusters = 100,
)
saveArchRProject(ArchRProj = proj, load = T)

png('~/yuzhao1/work/final_RC2atac/annotation/plots/stem/umap_loc.png',width = 1800, height = 2000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$biopsy_location
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/stem/umap_res0.5.png',width = 1800, height = 2000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res0.5
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/stem/umap_res1.png',width = 1800, height = 2000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res1
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/stem/umap_res1.5.png',width = 1800, height = 2000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$Harmony_Clusters_res1.5
plot_df_umap_custom(df, show.label = 'name')
dev.off()


proj$anno2 <- proj$anno1
proj$anno2[proj$Harmony_Clusters_res0.5 %in% c('C2', 'C3')] <- 'Stem1'
proj$anno2[proj$Harmony_Clusters_res0.5 %in% c('C1')] <- 'Stem2'

png('~/yuzhao1/work/final_RC2atac/annotation/plots/stem/umap_anno2.png',width = 1800, height = 2000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno2
plot_df_umap_custom(df, show.label = 'name')
dev.off()



saveArchRProject(ArchRProj = proj, load = T)









