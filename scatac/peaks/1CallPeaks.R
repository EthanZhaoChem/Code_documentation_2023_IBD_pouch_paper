library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(4)
our.dir <- '~/yuzhao1/work/manu/rc2/plots/4atac_MarkerPeaks/'

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj_tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj_bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj_myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj_others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")



############################ 1. epithelial ###############################

temp_lineage <- 'epithelial'  
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

# add pseudo group replicates
proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = "anno1",
                          useLabels = T,
                          minReplicates = 6,
                          maxReplicates = 24,
                          minCells = 50,
                          maxCells = 10000,
                          force = T)
## save
saveArchRProject(ArchRProj = proj, load = T)

# pseudo group peaks

pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "anno1",
  pathToMacs2 = pathToMacs2,
  reproducibility = "3",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 50,
  excludeChr = c("chrM", "chrY"),
  extsize = 150,
  cutOff = 0.05,
  extendSummits = 250,
  force = T
)
## save
saveArchRProject(ArchRProj = proj, load = T)

## add peak matrix to object
proj <- addPeakMatrix(proj, force = T)

## save
saveArchRProject(ArchRProj = proj, load = T)

getGroupBW(
  ArchRProj = proj,
  groupBy = "anno1",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE
)

proj$anno1_loc <- proj$anno1
proj$anno1_loc[proj$biopsy_location=='POU' & (proj$anno1 %in% c('EC1-1', 'EC1-2'))] <- 'EC-POU1'
proj$anno1_loc[proj$biopsy_location=='POU' & (proj$anno1 %in% c('EC2-1', 'EC2-2'))] <- 'EC-POU2'
proj$anno1_loc[proj$biopsy_location=='PP'  & (proj$anno1 %in% c('EC1-1', 'EC1-2', 'EC2-1', 'EC2-2'))] <- 'EC-PP'
proj$anno1_loc[proj$biopsy_location=='AC'  & (proj$anno1 %in% c('EC1-1', 'EC1-2', 'EC2-1', 'EC2-2'))] <- 'EC-AC'
proj$anno1_loc[proj$biopsy_location=='TI'  & (proj$anno1 %in% c('EC1-1', 'EC1-2', 'EC2-1', 'EC2-2'))] <- 'EC-TI'
saveArchRProject(ArchRProj = proj, load = T)

getGroupBW(
  ArchRProj = proj,
  groupBy = "anno1_loc",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE
)



############################ 2. tcell ###############################

temp_lineage <- 'tcell'  
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

# add pseudo group replicates
proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = "anno1",
                          useLabels = T,
                          minReplicates = 12,
                          maxReplicates = 24,
                          minCells = 40,
                          maxCells = 10000,
                          force = T)
## save
saveArchRProject(ArchRProj = proj, load = T)

# pseudo group peaks

pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "anno1",
  pathToMacs2 = pathToMacs2,
  reproducibility = "3",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 50,
  excludeChr = c("chrM", "chrY"),
  extsize = 150,
  cutOff = 0.05,
  extendSummits = 250,
  force = T
)
## save
saveArchRProject(ArchRProj = proj, load = T)

## add peak matrix to object
proj <- addPeakMatrix(proj, force = T)

## save
saveArchRProject(ArchRProj = proj, load = T)

# bigwig function can't identify +/-, use word
proj$anno1[proj$anno1=="CD103- CD4 Trm"] <- "CD103Negative CD4 Trm"
getGroupBW(
  ArchRProj = proj,
  groupBy = "anno1",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE
)





############################ 3. bcell ###############################

temp_lineage <- 'bcell'  
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

# add pseudo group replicates
proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = "anno1",
                          useLabels = T,
                          minReplicates = 6,
                          maxReplicates = 24,
                          minCells = 100,
                          maxCells = 10000,
                          force = T)
## save
saveArchRProject(ArchRProj = proj, load = T)

# pseudo group peaks

pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "anno1",
  pathToMacs2 = pathToMacs2,
  reproducibility = "3",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 50,
  excludeChr = c("chrM", "chrY"),
  extsize = 150,
  cutOff = 0.05,
  extendSummits = 250,
  force = T
)

## add peak matrix to object
proj <- addPeakMatrix(proj, force = T)

getGroupBW(
  ArchRProj = proj,
  groupBy = "anno1",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE
)

## save
saveArchRProject(ArchRProj = proj, load = T)





############################ 4. myeloid ###############################

temp_lineage <- 'myeloid'  
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

# add pseudo group replicates
proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = "anno1",
                          useLabels = T,
                          minReplicates = 6,
                          maxReplicates = 24,
                          minCells = 100,
                          maxCells = 10000,
                          force = T)

pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "anno1",
  pathToMacs2 = pathToMacs2,
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 50,
  excludeChr = c("chrM", "chrY"),
  extsize = 150,
  cutOff = 0.05,
  extendSummits = 250,
  force = T
)

## add peak matrix to object
proj <- addPeakMatrix(proj, force = T)


getGroupBW(
  ArchRProj = proj,
  groupBy = "anno1",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE
)

## save
saveArchRProject(ArchRProj = proj, load = T)



############################ 5. others ###############################

temp_lineage <- 'others'  
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)

# add pseudo group replicates
proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = "anno1",
                          useLabels = T,
                          minReplicates = 6,
                          maxReplicates = 24,
                          minCells = 60,
                          maxCells = 10000,
                          force = T)

pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "anno1",
  pathToMacs2 = pathToMacs2,
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 50,
  excludeChr = c("chrM", "chrY"),
  extsize = 150,
  cutOff = 0.05,
  extendSummits = 250,
  force = T
)

## add peak matrix to object
proj <- addPeakMatrix(proj, force = T)

getGroupBW(
  ArchRProj = proj,
  groupBy = "anno1",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE
)

## save
saveArchRProject(ArchRProj = proj, load = T)

