dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(harmony)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/figR/optMatching_functions.R')
source('~/yuzhao1/scripts/figR/utils.R')
library(optmatch)
library(Matrix)
library(FNN)
library(dplyr)
library(igraph)
library(pracma)
library(uwot) 
library(FigR)
library(foreach)

source('/project/gca/yuzhao1/work/manu/rc2/scripts/rna_annotation_markers.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')
source('~/yuzhao1/scripts/helper_archr.R')

proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
our.dir <- '~/yuzhao1/work/final_RC2atac/figr/epithelial_rds/'
######################### :) I am a divider :) ##########################
# #  get gene expression
# #  read files
# addArchRThreads(1)
# gim <- getMatrixFromProject(
#   ArchRProj = proj,
#   useMatrix = "GeneIntegrationMatrix",
#   useSeqnames = NULL,
#   verbose = TRUE,
#   binarize = FALSE,
#   threads = getArchRThreads(),
#   logFile = createLogFile("getMatrixFromProject")
# )
# 
# saveRDS(gim, paste0(proj@projectMetadata$outputDirectory, '/GeneIntegrationMatrix.rds'))


# the gene score is already scaled by total counts in each cell, should do log transform as the final normalization step
# the normalization functionin seurat will deal with it in the same way.
# gim <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/GeneIntegrationMatrix.rds'))
# 
# 
# counts <- gim@assays@data$GeneIntegrationMatrix
# rownames(counts) <- gim@elementMetadata$name
# colnames(counts) <- colnames(gim)
# idx <- match(colnames(gim), proj$cellNames)
# meta <- as.data.frame(proj@cellColData[idx,])
# 
# 
# seurat_rna <-
#   CreateSeuratObject(
#     counts,
#     project = "SeuratProject",
#     assay = "RNA",
#     min.cells = 0,
#     min.features = 0,
#     names.field = 1,
#     names.delim = "#",
#     meta.data = meta
#   )
# 
# DefaultAssay(seurat_rna) <- "RNA"
# seurat_rna <- NormalizeData(object = seurat_rna)
# 
# saveRDS(seurat_rna, paste0(proj@projectMetadata$outputDirectory, '/seurat_rna_normalized_GIM.rds'))


######################### :) I am a divider :) ##########################
addArchRThreads(1)
peak.mtx.all <- getMatrixFromProject(proj, useMatrix='PeakMatrix')
saveRDS(peak.mtx.all, paste0(proj@projectMetadata$outputDirectory, '/peak_matrix_unbinarized.rds'))

peak.mtx.all <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peak_matrix_unbinarized.rds'))
ATAC.se.paired <- peak.mtx.all[, proj$cellNames]
names(ATAC.se.paired@assays) <- 'counts'

# read rna files
seurat_rna <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/seurat_rna_normalized_GIM.rds'))
RNAmat.paired <- seurat_rna@assays$RNA@data[, proj$cellNames]


########################## :) I am a divider :) ##########################

# 1. run cisCorr
cisCorr <- runGenePeakcorr(ATAC.se = ATAC.se.paired,
                           RNAmat = RNAmat.paired,
                           genome = "hg38",
                           nCores = 24,
                           p.cut = NULL, # Set this to NULL and we can filter later
                           n_bg = 100)
saveRDS(cisCorr, paste0(our.dir, 'cisCorr.rds'))


########################## :) I am a divider :) ##########################

# 2. Filter peak-gene correlations by p-value
cisCor <- readRDS(paste0(our.dir, 'cisCorr.rds'))
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)


########################## :) I am a divider :) ##########################

# 3. Determine DORC genes (can't be run interactively)
dorcGenes <- cisCor.filt %>% dorcJPlot(cutoff=7, returnGeneList = TRUE)
saveRDS(dorcGenes, paste0(our.dir, 'dorcGenes.rds'))

# dorcGenes <- readRDS(paste0(our.dir, 'dorcGenes.rds'))
# human_tfs <- read.table('~/yuzhao1/resource/scenic/utoronto_human_tfs_v_1.01.txt', header = F)[[1]]
# intersect(dorcGenes, human_tfs)
# saveRDS(intersect(dorcGenes, human_tfs), '~/yuzhao1/work/final_RC2atac/figr/rds/dorcTFs.rds')

########################## :) I am a divider :) ##########################

# 4 Get DORC scores

peak.mtx.all <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peak_matrix_unbinarized.rds'))
ATAC.se.paired <- peak.mtx.all[, proj$cellNames]
names(ATAC.se.paired@assays) <- 'counts'

dorcMat <- getDORCScores(ATAC.se.paired, dorcTab=cisCor.filt, geneList=dorcGenes,nCores=24)
saveRDS(dorcMat, paste0(our.dir, 'dorcMat.rds'))

dorcMat <- readRDS(paste0(our.dir, 'dorcMat.rds'))

########################## :) I am a divider :) ##########################

# 5 TF-gene associations

# peak.mtx.all <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peak_matrix_unbinarized.rds'))
# ATAC.se.paired <- peak.mtx.all[, proj$cellNames]
# names(ATAC.se.paired@assays) <- 'counts'
# 
# # read rna files
# seurat_rna <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/seurat_rna_normalized_GIM.rds'))
# RNAmat.paired <- seurat_rna@assays$RNA@data[, proj$cellNames]
# 
# dorcMat <- dorcMat[, proj$cellNames]
# 
# 
# figR.d <- runFigRGRN(ATAC.se = ATAC.se.paired, # Must be the same input as used in runGenePeakcorr()
#                      dorcTab = cisCor.filt, # Filtered peak-gene associations
#                      genome = "hg38",
#                      dorcMat = dorcMat,
#                      rnaMat = RNAmat.paired,
#                      nCores = 8)
# saveRDS(figR.d, paste0(our.dir, 'figR.d.rds'))


# # visualization
# figR.d <- readRDS(paste0(our.dir, 'figR.d.rds'))
# 
# 
# # Visualize FigR results
# library(BuenColors)
# png(paste0(our.dir, '/DORC-to-TF associations.png'),res = 300, height = 1000, width = 1200)
# figR.d %>%
#   ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) +
#   ggrastr::geom_point_rast(size=0.01,shape=16) +
#   theme_classic() +
#   scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))
# dev.off()
# 
# rankDrivers(figR.d,rankBy = "meanScore",interactive = FALSE)
# 
# rankDrivers(figR.d,score.cut = 2,rankBy = "nTargets",interactive = TRUE)
# 
# plotDrivers(figR.d,score.cut = 1,marker = "SATB2")
# 
# 
# 
# 
# library(ComplexHeatmap)
# pdf(paste0(our.dir, 'figRHeatmap.pdf'), width = 10, height = 26)
# plotfigRHeatmap(figR.d = figR.d,
#                 score.cut = 1.5,
#                 column_names_gp = gpar(fontsize=12), # from ComplexHeatmap
#                 show_row_dend = FALSE # from ComplexHeatmap
# )
# dev.off()
# 
# # visually_filtered_ec_degs <- unique(c(rna_deg_markers_ec_ileum, rna_deg_markers_ec_colon))
# ec_degs <- readRDS('~/yuzhao1/work/manu/rc2/plots/2rna_deg/EC_heatmap/genesEC.rds')
# dorc.subset <- intersect(dorcGenes, ec_degs)
# pdf(paste0(our.dir, 'figRHeatmap_intersectDORC_ECdegRNA.pdf'), width = 10, height = 20)
# plotfigRHeatmap(figR.d = figR.d %>% subset(., DORC %in% dorc.subset),
#                 score.cut = 1.5,
#                 column_names_gp = gpar(fontsize=12), # from ComplexHeatmap
#                 show_row_dend = FALSE # from ComplexHeatmap
# )
# dev.off()
# 
# 
# 
# 
# library(networkD3)
# png(paste0(our.dir, 'figR_network.png'),res = 300, height = 3000, width = 3000)
# plotfigRNetwork(figR.d = figR.d,
#                 score.cut = 2,
#                 TFs = NULL,
#                 weight.edges = TRUE)
# dev.off()
# 
# 
# png(paste0(our.dir, 'figR_network_intersectDORC_ECdegRNA.png'),res = 300, height = 6000, width = 3000)
# plotfigRNetwork(figR.d = figR.d %>% subset(., DORC %in% dorc.subset),
#                 score.cut = 2,
#                 TFs = NULL,
#                 weight.edges = TRUE)
# dev.off()
# 
# 
# 
# # # ########################### :) I am a divider :) ##########################
# # 11 add dorc assay to the paired rna seurat subset
# seurat_rna <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/seurat_rna_normalized_GIM.rds'))
# RNAmat.paired <- seurat_rna@assays$RNA@data[, proj$cellNames]
# dorcMat <- dorcMat[, proj$cellNames]
# 
# # add updated centered dorc counts to seurat "DORC" assay
# seurat_rna[["DORC"]] <- CreateAssayObject(data = dorcMat)
# DefaultAssay(seurat_rna) <- "DORC"
# 
# 
# # # ########################### :) I am a divider :) ##########################
# # 12. Find differential testing of DORC accessibility scores or expression levels, use top 10 from DE DORC score
# # it means more regions in one group can regulate this TF, this is all what it means
# 
# EC.markers1 <- FindMarkers(seurat_rna, ident.1 = 'EC-AC', ident.2 = c("EC-TI"), group.by = 'anno1.loc', logfc.threshold=1)
# EC.markers2 <- FindMarkers(seurat_rna, ident.1 = 'EC-POU2', ident.2 = c("EC-PP"), group.by = 'anno1.loc', logfc.threshold=1)
# 
# 
# VlnPlot(seurat_rna, features = 'SATB2', pt.size = 0, group.by = 'anno1.loc')
# VlnPlot(seurat_rna, features = 'TFCP2L1', pt.size = 0, group.by = 'anno1.loc')
# VlnPlot(seurat_rna, features = 'MAF', pt.size = 0, group.by = 'anno1.loc')
# VlnPlot(seurat_rna, features = 'CA2', pt.size = 0, group.by = 'anno1.loc')
# VlnPlot(seurat_rna, features = 'SATB2', pt.size = 0, group.by = 'anno1.loc')
# 
# # # ########################### :) get motif pfm :) ##########################
# 
# packagePath <- find.package("FigR", lib.loc=NULL, quiet = TRUE)
# 
# pwm <- readRDS(paste0(packagePath,"/data/cisBP_human_pfms_2021.rds"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# xx <- read.table('~/yuzhao1/work/final_GCArna/metadata/meta_Ethan_curated_20221108.csv', header =T, sep = ',')
# 
# yy <- list.files('~/gca/GCA_scRNA/fastq')
# 
# 
# unique(xx$Patient_ID)





















