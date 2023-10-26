dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(4)

# # epithelial
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  groupBy = "anno2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno2_epithelial', '.rds'))



# # tcell
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  groupBy = "anno1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))


# # bcell
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  groupBy = "anno1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))


# # myeloid
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  groupBy = "anno1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))

# # stromal
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  groupBy = "anno1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1', '.rds'))




