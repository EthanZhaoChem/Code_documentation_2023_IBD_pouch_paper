library(Seurat)
library(arrow)
library(progressr)
library(udunits2)
library(sf)
source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')

n_samples <- 6
paths_samples <- c('ac' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__50AC__20240711__171029',
                   'ti' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__50TI__20240711__171029',
                   'pouR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__72pov-222__20240711__171029',
                   'ppR1' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036222__72pp-222__20240711__171029',
                   'pouR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__72pov-314__20240711__171029',
                   'ppR2' = '/project/spott/collaborations/Pouch_Xenium/Run_1/Output/output-XETG00120__0036314__72pp-314__20240711__171029')

xenium.obj.list <- list()
for (i in 1:n_samples) {
  sample_name <- names(paths_samples)[[i]]
  tmp.obj <- LoadXenium_alt(paths_samples[[sample_name]], fov = "fov")
  xenium.obj.list[[sample_name]] <- subset(tmp.obj, subset = nCount_Xenium > 0)
}

##################  ################ ##################  ################
# # QC
# VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

##################  ################ ##################  ################
# normalization, dimensionality reduction
xenium.obj <- xenium.obj.list[['ac']]
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.1)

DimPlot(xenium.obj)
FeaturePlot(xenium.obj, features = c("CEACAM5", "CDX2", "FABP6"))
ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)


##################  ################ ##################  ################
# plot whole image
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("CEACAM5", "CDX2", "FABP6"), nmols = 20000)
ImageFeaturePlot(xenium.obj,
                 features = c("CEACAM5", "CDX2", "FABP6"),
                 min.cutoff = c(0, 0, 0),
                 max.cutoff = c(1,5,30),
                 size = 0.75, cols = c("white", "red"))

##################  ################ ##################  ################
# plot a ROI
cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 2000), y = c(0, 3000), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white",
             border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,
             molecules = c("CEACAM5", "CDX2", "FABP6"), nmols = 10000)
























