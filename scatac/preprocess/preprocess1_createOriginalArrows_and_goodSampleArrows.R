library(stringr)
library(ArchR)
library(Seurat)
# source('~/yuzhao1/scripts/plot.R')

# read files
setwd('~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2/')
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2/projdir/")


# set global parameters
filterRatio = 2.0
minFrags <- 2000
minTSS <- 6
GenomeSet <- 'hg38'
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
addArchRGenome(GenomeSet)
addArchRThreads(12)
inputFiles <- c('HA01-AC_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-HERZ-8S-scATAC-HA01-AC_deep/fragments.tsv.gz",
                'HA01-TI_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-HERZ-8S-scATAC-HA01-TI_deep/fragments.tsv.gz",
                'HA02-AC_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-HERZ-8S-scATAC-HA02-AC_deep/fragments.tsv.gz", 
                'HA02-TI_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-HERZ-8S-scATAC-HA02-TI_deep/fragments.tsv.gz", 
                'HA04-AC_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-HERZ-8S-scATAC-HA04-AC_deep/fragments.tsv.gz", 
                'HA04-TI_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-HERZ-8S-scATAC-HA04-TI_deep/fragments.tsv.gz", 
                'HA50-TI_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-HA50-TI_deep/fragments.tsv.gz",
                'HA50-AC_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AAB-RZ-RZ20210809-HA50-AC_deep/fragments.tsv.gz", 
                'HA51-AC_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-HA51-AC_deep/fragments.tsv.gz", 
                'HA51-TI_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-HA51-TI_deep/fragments.tsv.gz", 
                'HA55-AC_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-HA55-AC_deep/fragments.tsv.gz", 
                'HA55-TI_deep' = "~/gca/GCA_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-HA55-TI_deep/fragments.tsv.gz", 
                'OR72-POU' = "/project/gca/Pouch_scATAC/data/CellRanger_Output/AB-RZ-10X-12S-ATAC-20220111S1-OR72POU/fragments.tsv.gz", 
                'OR72-PP' = "/project/gca/Pouch_scATAC/data/CellRanger_Output/AB-RZ-10X-12S-ATAC-20220111S1-OR72PP/fragments.tsv.gz", 
                'OR43-POU' = "/project/gca/Pouch_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-OR43-POU/fragments.tsv.gz", 
                'OR43-PP' = "/project/gca/Pouch_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-OR43-PP/fragments.tsv.gz", 
                'OR48-POU' = "/project/gca/Pouch_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-OR48-POU/fragments.tsv.gz", 
                'OR48-PP' = "/project/gca/Pouch_scATAC/data/CellRanger_Output/AB-RZ-RZ20210809-OR48-PP/fragments.tsv.gz", 
                'OR101-POU' = "/project/gca/Pouch_scATAC/data/CellRanger_Output20220417/AB-RZ-0321-OR101-POU/outs/fragments.tsv.gz", 
                'OR101-PP' = "/project/gca/Pouch_scATAC/data/CellRanger_Output20220417/AB-RZ-0321-OR101-PP/outs/fragments.tsv.gz", 
                'OR102-POU' = "/project/gca/Pouch_scATAC/data/CellRanger_Output20220417/AB-RZ-0321-OR102-POU/outs/fragments.tsv.gz", 
                'OR102-PP' = "/project/gca/Pouch_scATAC/data/CellRanger_Output20220417/AB-RZ-0321-OR102-PP/outs/fragments.tsv.gz", 
                'OR109-POU' = "/project/gca/Pouch_scATAC/data/CellRanger_Output20220417/AB-RZ-0321-OR109-POU/outs/fragments.tsv.gz", 
                'OR109-PP' = "/project/gca/Pouch_scATAC/data/CellRanger_Output20220417/AB-RZ-0321-OR109-PP/outs/fragments.tsv.gz")







#prepare arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  outputNames = names(inputFiles),
  minTSS = minTSS,
  minFrags = minFrags,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  verbose = T,
  force = T
)


# use arrowfiles directly
ArrowFiles <- paste0(names(inputFiles), '.arrow')

# # infer doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  nTrials = 10,
  dimsToUse = 2:50,
  force = T
)

# create, clean and reload archr project
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "projdir", copyArrows = F)
proj <- filterDoublets(ArchRProj = proj, filterRatio = filterRatio)
saveArchRProject(ArchRProj = proj, load = T)
















