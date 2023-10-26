dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(4)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'

proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
# proj <- addImputeWeights(proj)
# saveArchRProject(ArchRProj = proj, load = T)

# ############################# Section extra: anno2 update ################################
proj_stem <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_Stem/")
cells1 <- proj_stem$cellNames[proj_stem$anno2 == 'Stem1'] 
cells2 <- proj_stem$cellNames[proj_stem$anno2 == 'Stem2'] 

proj$anno2 <- proj$anno1
idx1 <- match(cells1, proj$cellNames)
proj$anno2[idx1] <- 'Stem1' 
idx2 <- match(cells2, proj$cellNames)
proj$anno2[idx2] <- 'Stem2' 

saveArchRProject(ArchRProj = proj, load = T)


# ############################# Section extra: anno1.loc update ################################
proj$anno1.loc <- paste0(proj$anno1, '-', proj$biopsy_location)
proj$anno1.loc[proj$anno2 %in% c("Stem1") & proj$biopsy_location =='POU'] <- 'Stem-POU1'
proj$anno1.loc[proj$anno2 %in% c("Stem2") & proj$biopsy_location =='POU'] <- 'Stem-POU2'

proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2") & proj$biopsy_location =='POU'] <- 'EC-POU1'
proj$anno1.loc[proj$anno1 %in% c("EC2-1", "EC2-2") & proj$biopsy_location =='POU'] <- 'EC-POU2'
proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2") & proj$biopsy_location =='PP'] <- 'EC-PP'
proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2") & proj$biopsy_location =='TI'] <- 'EC-TI'
proj$anno1.loc[proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2") & proj$biopsy_location =='AC'] <- 'EC-AC'

proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='POU'] <- 'Goblet-POU'
proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='PP'] <- 'Goblet-PP'
proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='TI'] <- 'Goblet-TI'
proj$anno1.loc[proj$anno1 %in% c("Goblet1", "Goblet2") & proj$biopsy_location =='AC'] <- 'Goblet-AC'

proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='POU'] <- 'CD4T-POU'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='PP'] <- 'CD4T-PP'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='TI'] <- 'CD4T-TI'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD4 Trm", "CD103- CD4 Trm", "CD4 Tcm", "Treg") & proj$biopsy_location =='AC'] <- 'CD4T-AC'

proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='POU'] <- 'CD8T-POU'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='PP'] <- 'CD8T-PP'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='TI'] <- 'CD8T-TI'
proj$anno1.loc[proj$anno1 %in% c("CD103+ CD8 Trm", "KLRG1+ CD8 Trm") & proj$biopsy_location =='AC'] <- 'CD8T-AC'

proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='POU'] <- 'Stromal-POU'
proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='PP'] <- 'Stromal-PP'
proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='TI'] <- 'Stromal-TI'
proj$anno1.loc[proj$anno1 %in% c("Stromal-1", "Stromal-2", "Stromal-3") & proj$biopsy_location =='AC'] <- 'Stromal-AC'

proj$anno1_loc <- proj$anno1.loc
saveArchRProject(ArchRProj = proj, load = T)

# ############################# Section extra: anno1_loc update ################################

proj_union <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj_immune <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_immune_filtered/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj_tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj_bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj_myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj_others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")

proj_union$anno1_loc <- proj_union$anno1.loc
proj_immune$anno1_loc <- proj_immune$anno1.loc
proj_epithelial$anno1_loc <- proj_epithelial$anno1.loc
proj_tcell$anno1_loc <- proj_tcell$anno1.loc
proj_bcell$anno1_loc <- proj_bcell$anno1.loc
proj_myeloid$anno1_loc <- proj_myeloid$anno1.loc
proj_others$anno1_loc <- proj_others$anno1.loc

saveArchRProject(ArchRProj = proj_union, load = T)
saveArchRProject(ArchRProj = proj_immune, load = T)
saveArchRProject(ArchRProj = proj_epithelial, load = T)
saveArchRProject(ArchRProj = proj_tcell, load = T)
saveArchRProject(ArchRProj = proj_bcell, load = T)
saveArchRProject(ArchRProj = proj_myeloid, load = T)
saveArchRProject(ArchRProj = proj_others, load = T)




# ############################### Section0: make filtered union set ################################
proj.epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj.tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj.bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj.myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj.others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")
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
cellsPass <- c(proj.epithelial$cellNames, proj.tcell$cellNames, proj.bcell$cellNames, proj.myeloid$cellNames, proj.others$cellNames)
cellsAnno1 <- c(proj.epithelial$anno1, proj.tcell$anno1, proj.bcell$anno1, proj.myeloid$anno1, proj.others$anno1)
proj@cellColData$anno1 <- mapvalues(proj$cellNames, from = cellsPass, to = cellsAnno1)%>% unlist()
saveArchRProject(ArchRProj = proj, load = T)

# add predicted labels to the filtered proj
cellsPass <- c(proj.epithelial$cellNames, proj.tcell$cellNames, proj.bcell$cellNames, proj.myeloid$cellNames, proj.others$cellNames)
labels_predicted <- c(proj.epithelial$predictedGroup_lineage_constrained, proj.tcell$predictedGroup_lineage_constrained, 
                      proj.bcell$predictedGroup_lineage_constrained, proj.myeloid$predictedGroup_lineage_constrained,
                      proj.others$predictedGroup_lineage_constrained)
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
  varFeatures = 150000,
  clusterParams = list(
    sampleCells = 60000,
    resolution = c(0.5, 1, 1.5, 2, 2),
    n.start = 10,
    maxClusters = 15
  ),
  sampleCellsPre = 60000,
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


png('~/yuzhao1/work/final_RC2atac/annotation/plots/union_filtered/umap_anno1.png',width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png('~/yuzhao1/work/final_RC2atac/annotation/plots/union_filtered/umap_anno1_loc.png',width = 5000, height = 5200,res = 300)
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






