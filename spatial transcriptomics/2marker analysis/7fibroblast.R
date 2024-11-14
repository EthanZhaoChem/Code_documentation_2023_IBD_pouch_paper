library(Seurat)
library(arrow)
library(progressr)
library(udunits2)
library(sf)
library(patchwork)
source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')

# read xenium
n_samples <- 8
paths_samples <- c('ac' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output//output-XETG00120__0036222__50AC__20240711__171029',
                   'ti' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output//output-XETG00120__0036222__50TI__20240711__171029',
                   'pouR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output//output-XETG00120__0036222__72pov-222__20240711__171029',
                   'ppR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output//output-XETG00120__0036222__72pp-222__20240711__171029',
                   'pouR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output//output-XETG00120__0036314__72pov-314__20240711__171029',
                   'ppR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output//output-XETG00120__0036314__72pp-314__20240711__171029',
                   'pou53' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__53POV__20240711__171029',
                   'pp53' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__53PP__20240711__171029')

xenium.obj.list <- list()
for (i in 1:n_samples) {
  sample_name <- names(paths_samples)[[i]]
  tmp.obj <- LoadXenium_alt(paths_samples[[sample_name]], fov = "fov")
  xenium.obj.list[[sample_name]] <- subset(tmp.obj, subset = nCount_Xenium > 0)
}

# read seurat
seurat <- readRDS('~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged_processed.rds')
metadata <- seurat@meta.data

# set path
out.dir <- '~/yuzhao1/work/final_RC2rna/0revision/spatial2/plots/fibroblast/'

##################  ################ ##################  ################
# pre select a few genes to plot
genes_selected <- c('CEACAM5', 'FABP6', 'ADAMDEC1', 'NRG1', 'GREM1', 'OGN')
max_cutoffs <- rep(3, 5)
min_cutoffs <- c(0, 0, 0, 0, 0)
p2.width <- 50

##################  ################ ##################  ################
# 1. pouch
tmp.sample <- 'pouR1'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(400, 1000), x = c(1400, 2000), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom', 
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()



##################  ################ ##################  ################
# 2. pp
tmp.sample <- 'ppR1'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(600, 1200), x = c(600, 1200), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom',
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()


##################  ################ ##################  ################
# 3. ti
tmp.sample <- 'ti'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(400, 1600), x = c(500, 1400), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom',
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()



##################  ################ ##################  ################
# 4. ac
tmp.sample <- 'ac'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(1500, 2100), x = c(0, 600), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom',
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()



##################  ################ ##################  ################
# 5. pou53 
tmp.sample <- 'pou53'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(500, 1100), x = c(1300, 1900), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom',
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()



##################  ################ ##################  ################
# 6. pp53 
tmp.sample <- 'pp53'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(300, 900), x = c(800,1400), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom',
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()



##################  ################ ##################  ################
# 7. pouch 72 replicate 2 
tmp.sample <- 'pouR2'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(600,1400), x = c(1800,2500), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom',
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()

##################  ################ ##################  ################
# 8. pp 72 replicate 2 
tmp.sample <- 'ppR2'
xenium.obj <- xenium.obj.list[[tmp.sample]]
xenium.obj$predicted.celltype <- metadata[paste0(tmp.sample, '_', Cells(xenium.obj)), "predicted.celltype"]
cropped.coords <- Crop(xenium.obj[["fov"]], y = c(250, 850), x = c(1300,1900), coords = "plot")

xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
p1 <- ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
                   border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
                   molecules = genes_selected, nmols = 10000)+
  plot_layout(nrow = 1)

p2 <- ImageFeaturePlot(xenium.obj,fov = 'zoom',
                       features = genes_selected,
                       min.cutoff = min_cutoffs,
                       max.cutoff = max_cutoffs,
                       size = 0.75, 
                       cols = c("white", "red"))+ 
  plot_layout(nrow = 1)

pdf(paste0(out.dir, tmp.sample, '_part1.pdf'), height = 5, width = 6)
print(p1)
dev.off()

pdf(paste0(out.dir, tmp.sample, '_part2.pdf'), height = 8, width = p2.width)
print(p2)
dev.off()

