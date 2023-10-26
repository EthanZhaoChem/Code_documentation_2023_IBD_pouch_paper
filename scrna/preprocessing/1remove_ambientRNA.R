dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(Seurat)
library(SoupX)
library(Matrix)
library(foreach)
setwd('~/gca/yuzhao1/work/final_RC2rna/upstream')
SampleIDs <- read.csv("~/yuzhao1/work/final_RC2rna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID
SampleIDs <- c(SampleIDs, 'OR53-POU', 'OR53-PP')
SampleIDs <- SampleIDs[c(grep('POU', SampleIDs), grep('PP', SampleIDs))]


mtx_outdir_removedAmbientRNA <- '~/yuzhao1/work/final_RC2rna/preprocessing/matrix_removedAmbientRNA/'
gene_outdir_removedAmbientRNA <- '~/yuzhao1/work/final_RC2rna/preprocessing/gene_removedAmbientRNA/'
barcodes_outdir <- '~/yuzhao1/work/final_RC2rna/preprocessing/barcodes_outdir/'


# backend (need to customize max iterations)
max_iterations <- 14
opts <- list()
pb <- txtProgressBar(min = 0, max = max_iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# set parallel backend
nCores = 7
cl <- parallel::makeCluster(nCores)
doSNOW::registerDoSNOW(cl)
getDoParWorkers()
getDoParName()

# do parallel computing
require(foreach)
foreach::foreach(i=1:max_iterations , .packages = c('patchwork', "dplyr", "Matrix", "SoupX", "Seurat"), .options.snow = opts, .inorder = F) %dopar% {
  # no need to re-process AC-TI data
  if(grepl('TI', SampleIDs[[i]]) | grepl('AC', SampleIDs[[i]])){
    next
  }
  # Load data and estimate soup profile
  cellranger_outFolder <- paste0(SampleIDs[[i]], '/outs')
  sc = load10X(cellranger_outFolder)
  # Estimate rho
  sc = autoEstCont(sc)
  # Clean the data
  out = adjustCounts(sc)
  # Add SampleID to cell barcodes to make them unique
  colnames(out) = paste0(SampleIDs[[i]], "_", colnames(out))
  # Store the result
  writeMM(out, paste0(mtx_outdir_removedAmbientRNA, SampleIDs[[i]], '.mtx'))
  write(x = rownames(out), file = paste0(gene_outdir_removedAmbientRNA , SampleIDs[[i]], '.tsv'))
  write(x = colnames(out), file = paste0(barcodes_outdir, SampleIDs[[i]], '.tsv'))
  
  srat[[SampleIDs[[i]]]] = out
}
parallel::stopCluster(cl)



