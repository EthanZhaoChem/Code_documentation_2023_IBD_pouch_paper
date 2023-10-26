dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(4)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'

proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_immune_filtered/")
# proj <- addImputeWeights(proj)
# saveArchRProject(ArchRProj = proj, load = T)





# ############################# Section extra: anno1.loc update ################################
proj$anno1.loc <- paste0(proj$anno1, '-', proj$biopsy_location)

proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='POU'] <- 'CD4T-POU'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='PP'] <- 'CD4T-PP'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='TI'] <- 'CD4T-TI'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='AC'] <- 'CD4T-AC'

proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='POU'] <- 'CD8T-POU'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='PP'] <- 'CD8T-PP'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='TI'] <- 'CD8T-TI'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='AC'] <- 'CD8T-AC'

proj$anno1_loc <- proj$anno1.loc
saveArchRProject(ArchRProj = proj, load = T)


# ############################### Section0: make filtered union set ################################
proj.tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj.bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj.myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
# 
# cellsPass <- c(proj.epithelial$cellNames, proj.tcell$cellNames, proj.bcell$cellNames, proj.myeloid$cellNames, proj.others$cellNames)
# proj.filtered <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = cellsPass,
#   outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered",
#   dropCells = T,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = T
# )
# 
# cellsPass.immune <- c(proj.tcell$cellNames, proj.bcell$cellNames, proj.myeloid$cellNames)
# proj.filtered.immune <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = cellsPass.immune,
#   outputDirectory = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_immune_filtered",
#   dropCells = TRUE,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = T
# )

# add anno1 labels to the filtered proj
cellsPass <- c(proj.tcell$cellNames, proj.bcell$cellNames, proj.myeloid$cellNames)
cellsAnno1 <- c(proj.tcell$anno1, proj.bcell$anno1, proj.myeloid$anno1)
proj@cellColData$anno1 <- mapvalues(proj$cellNames, from = cellsPass, to = cellsAnno1)%>% unlist()
saveArchRProject(ArchRProj = proj, load = T)

# add predicted labels to the filtered proj
cellsPass <- c(proj.tcell$cellNames, proj.bcell$cellNames, proj.myeloid$cellNames)
labels_predicted <- c(proj.tcell$predictedGroup_lineage_constrained, 
                      proj.bcell$predictedGroup_lineage_constrained, 
                      proj.myeloid$predictedGroup_lineage_constrained)

proj@cellColData$predictedGroup_lineage_constrained <- mapvalues(proj$cellNames, from = cellsPass, to = labels_predicted)%>% unlist()
saveArchRProject(ArchRProj = proj, load = T)




# ############################# Section1: Dim-Red ################################
# Dim Reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 6,
  dimsToUse = 2:100,
  varFeatures = 100000,
  clusterParams = list(
    sampleCells = 20000,
    resolution = c(0.5, 1, 1.5, 2, 2),
    n.start = 10,
    maxClusters = 10
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
  minDist = 0.2,
  force = T,
  dimsToUse = 1:99,
  metric = "cosine"
)

saveArchRProject(ArchRProj = proj, load = T)


png('~/yuzhao1/work/final_RC2atac/annotation/plots/union_filtered_immune/umap_anno1.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/union_filtered_immune/umap_anno1_loc.png',width = 5000, height = 5200,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'name')+
    facet_wrap(~ biopsy_location) +
    theme(
      strip.background = element_rect(fill = "white", colour = "white"),
      strip.text = element_text(size = 12)
    )
dev.off()






