dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(4)
# our.dir <- '~/yuzhao1/work/manu/rc2/plots/4atac_MarkerPeaks/'

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_filtered/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
proj_tcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_tcell_filtered2/")
proj_bcell <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_bcell_filtered/")
proj_myeloid <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_myeloid_filtered2/")
proj_others <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_others_filtered2/")

##################### prepare database ############################
# # download file from JASPAR2022
# library(TFBSTools)
# library(JASPAR2022)
# # read TF databases
# pwm_set1 <- getMatrixSet(x = JASPAR2022, opts = list(all_versions = FALSE, species = 'Homo sapiens', collection="CORE", matrixtype="PWM"))
# pwm_set2 <- getMatrixSet(x = JASPAR2022, opts = list(all_versions = FALSE, species = 'Homo sapiens', collection="UNVALIDATED", matrixtype="PWM"))
# 
# # give a unique TF name to avoid duplicates from unvalidated database
# names(pwm_set1) <- name(pwm_set1)
# names(pwm_set2) <- paste0(name(pwm_set2), '_unvalidated')
# pwm_set <- c(pwm_set1, pwm_set2)
# saveRDS(pwm_set, '~/yuzhao1/resource/JASPAR2022_Homosapiens_Core_and_unvalidated.rds')
# 
# pwm_set <- readRDS('~/yuzhao1/resource/JASPAR2022_Homosapiens_Core_and_unvalidated.rds')
# names(pwm_set) <- make.unique(names(pwm_set), sep = "._.") #adds a unique suffix onto end of duplicates


############################ 1. calculation ###############################
# for (temp_lineage in unique(proj_union$lineage)) {
#   proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
# #   proj <- addMotifAnnotations(ArchRProj = proj, motifPWMs = pwm_set, force = TRUE)
#   proj <- addMotifAnnotations(ArchRProj = proj,  motifSet = "cisbp", force = TRUE)
#   proj <- addBgdPeaks(proj, force = TRUE)
#   ## save
#   saveArchRProject(ArchRProj = proj, load = T)
# }
# 
# 
# for (temp_lineage in unique(proj_union$lineage)) {
#   proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
#   proj <- addDeviationsMatrix(
#     ArchRProj = proj,
#     # peakAnnotation = 'Motif',
#     matrixName = 'cisbp',
#     out = c("z", "deviations"),
#     binarize = FALSE,
#     verbose = TRUE,
#     force = T
#   )
#   saveArchRProject(ArchRProj = proj, load = T)
#   
#   ## marker peaks
#   markersPeaks <- getMarkerFeatures(
#     ArchRProj = proj, 
#     useMatrix = "PeakMatrix", 
#     normBy = 'ReadsInTSS',
#     maxCells = 20000,
#     groupBy = "anno1",
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     testMethod = "wilcoxon"
#   )
#   
#   saveRDS(markersPeaks, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_', temp_lineage, '.rds'))
# 
# }
# 
# for (temp_lineage in unique(proj_union$lineage)) {
#   proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
#   
#   markersPeaks <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_', temp_lineage, '.rds'))
#   enrichMotifs <- peakAnnoEnrichment(
#     seMarker = markersPeaks,
#     ArchRProj = proj,
#     peakAnnotation = "Motif",
#     cutOff = "FDR <= 0.01 & Log2FC >= 1"
#   )
#   
#   saveRDS(enrichMotifs, paste0(proj@projectMetadata$outputDirectory, '/', 'enrichMotifs_', temp_lineage, '.rds'))
# }

# group EC and goblet in epithelial
temp_lineage <- 'epithelial'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj$anno1.grouped <- proj$anno1
proj$anno1.grouped[proj$anno1 %in% c("EC1-1", "EC1-2", "EC2-1", "EC2-2")] <- 'EC'
proj$anno1.grouped[proj$anno1 %in% c("Goblet1", "Goblet2")] <- 'Goblet'

markersPeaks_grouped <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  normBy = 'ReadsInTSS',
  maxCells = 20000,
  groupBy = "anno1.grouped",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks_grouped, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_grouped_', temp_lineage, '.rds'))

enrichMotifs_grouped <- peakAnnoEnrichment(
  seMarker = markersPeaks_grouped,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

saveRDS(enrichMotifs_grouped, paste0(proj@projectMetadata$outputDirectory, '/', 'enrichMotifs_grouped_', temp_lineage, '.rds'))



################# save imputated cisbp matrix #############
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")


cisbp.mtx <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "cisbp",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

cisbp.mtx.z <- cisbp.mtx@assays@data$z
colnames(cisbp.mtx.z) <- rownames(cisbp.mtx@colData)

cisbp.mtx.z.imputed <- imputeMatrix(
  mat = cisbp.mtx.z,
  imputeWeights = getImputeWeights(proj),
  threads = 1,
  verbose = FALSE,
  logFile = createLogFile("imputeMatrix")
)

saveRDS(cisbp.mtx.z.imputed, '~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/ChromVAR_cisbp.mtx.z.imputed.rds')

