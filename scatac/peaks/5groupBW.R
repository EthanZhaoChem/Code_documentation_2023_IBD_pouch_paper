library(stringr)
library(ArchR)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(1)

# # read files
# proj_union <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
# proj_tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
# proj_bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
# proj_myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
# proj_others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")


############################ 1. union ###############################

temp_lineage <- 'epithelial'  
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)



getGroupBW(
  ArchRProj = proj,
  groupBy = "anno2",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 100000,
  ceiling = 4,
  verbose = TRUE
)

