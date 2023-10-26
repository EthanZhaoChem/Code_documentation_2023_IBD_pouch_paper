dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(limma)
library(edgeR)
library(variancePartition)
library(BiocParallel)
library(tidyverse)

source('~/yuzhao1/scripts/helper_archr.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')
addArchRThreads(4)
# our.dir <- '~/yuzhao1/work/manu/rc2/plots/4atac_MarkerPeaks/'

# # read files
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
metadata.original <- as.data.frame(proj@cellColData)

peakMtx.s4 <-  readRDS('~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/peak_matrix_unbinarized.rds')
peakMtx <- peakMtx.s4@assays@data$PeakMatrix
xx <- peakMtx.s4@rowRanges
rownames(peakMtx) <- paste0(seqnames(xx), '_', start(xx), '_', end(xx))

### POU_EC2vsPOU_EC1
metadata.original$anno1[metadata.original$anno1 %in% c("EC1-1", "EC1-2")] <- 'EC1'
metadata.original$anno1[metadata.original$anno1 %in% c("EC2-1", "EC2-2")] <- 'EC2'
metadata.original$anno1.sub <- paste0(metadata.original$biopsy_location,'_',metadata.original$anno1)
metadata.original <- subset(metadata.original, anno1 == 'EC1'| anno1 == 'EC2')
metadata.original$anno1[metadata.original$anno1 %in% c("EC1", "EC2")] <- 'EC'
unique(metadata.original$anno1)
unique(metadata.original$anno1.sub)
metadata.original$anno1.sub[metadata.original$anno1.sub %in% c("PP_EC1", "PP_EC2")] <- 'PP_EC'
metadata.original <- metadata.original[metadata.original$anno1.sub %in% c("PP_EC", "POU_EC2"),]

# design formula
classification = 'anno1'
model_formula = ~  anno1.sub + (1|patient)
maineffect = 'anno1.sub'
pseudo_factors = c('anno1.sub', 'patient')
cells.use <- rownames(metadata.original)

# input
object <- peakMtx[,cells.use]
labels <- pseudo_factors

factorlist <- list()
for (i in labels) factorlist[[i]] <- unique(metadata.original[,i])

meta <- expand.grid(factorlist, stringsAsFactors = FALSE)
rownames(meta) <- apply(meta, 1, function(x) paste0(x, collapse = '.'))

n <- nrow(meta)
out <- matrix(nrow=dim(object)[1], ncol=n, data=0)
rownames(out) <- rownames(object)
colnames(out) <- rownames(meta)

ncells <- c()
ncounts <- c()
total.cells <- dim(object)[2]
for (i in 1:n){
  cells <- 1:total.cells
  for (j in names(meta)) {
    keep  <- which(metadata.original[[j]] == meta[i,j])
    cells <- cells[cells %in% keep]
  }
  ncells[i] <- length(cells)
  ncounts[i] <- sum(object[,cells])
  
  #some other thing to measure
  if (length(cells)==1) {
    out[,i] <- object[,cells]
  } else {
    out[,i] <- Matrix::rowSums(object[,cells]) # all counts for that gene in this pesudo bulk
  }
}
meta$ncells <- ncells
meta$ncounts <- ncounts/max(ncounts) # normalize by max counts
#add that something else as metadata
pseudo <- list(counts=out, meta=meta)


# filter bulks
threshold <- 200
w <- which(pseudo$meta$ncells > threshold)
pseudo$counts <- pseudo$counts[,w]
pseudo$meta <- pseudo$meta[w,]

# differential test: should be at least 3 pesudo bulks
if (ncol(pseudo$counts) < 3){
  next
}

d <- DGEList(pseudo$counts)
d$samples <- cbind(d$samples[,c("lib.size","norm.factors")], pseudo$meta)

# # don't filter peaks because it is already filtered in archr
# keepgenes <- filterByExpr(d$counts, group = d$samples[[maineffect]])
# d <- d[keepgenes,]

d <- calcNormFactors(d, method = "TMM")
v <- voomWithDreamWeights(d, model_formula, d$samples, plot=FALSE)
modelfit <- dream(exprObj = v, formula = model_formula, data = d$samples, quiet = F, suppressWarnings = TRUE)
results <- modelfit
saveRDS(results ,'~/yuzhao1/work/final_RC2atac/dar_great/rds/da_results_POU_EC2vsPP.rds')


# # 7. analyze dar results
# de_results <- readRDS('~/yuzhao1/work/final_RC2atac/dar_great/rds/da_results_AC_ECvsTI_EC.rds')
# 
# fclist <- topTable(
#   de_results[["EC"]],
#   coef = "biopsy_locationTI",
#   p.value = 0.001,
#   sort.by = 'logFC',
#   number = 30000
# )





