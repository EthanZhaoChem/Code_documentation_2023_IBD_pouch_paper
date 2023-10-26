dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(harmony)
library(Seurat)
library(limma)
library(edgeR)
library(variancePartition)
library(BiocParallel)
library(tidyverse)

source('~/yuzhao1/scripts/plot.R')
source('/project/gca/yuzhao1/work/manu/rc2/scripts/rna_annotation_markers.R')
source('/project/gca/yuzhao1/work/manu/rc2/scripts/umap_colors.R')
source('/project/gca/yuzhao1/scripts/seurat/deg_pseudobulk.R')
out.dir <- '/project/gca/yuzhao1/work/manu/rc2/plots/'


epithelial <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/epithelial3.rds')
epithelial <- subset(epithelial, biopsy_location == 'POU')


### POU_EC2vsPOU_EC1
seurat <- epithelial
unique(seurat$anno1)

seurat$anno1[seurat$anno1 %in% c("EC1-1", "EC1-2")] <- 'EC1'
seurat$anno1[seurat$anno1 %in% c("EC2-1", "EC2-2")] <- 'EC2'
seurat$anno1.sub <- paste0(seurat$biopsy_location,'_',seurat$anno1)
seurat <- subset(seurat, anno1 == 'EC1'| anno1 == 'EC2')
seurat$anno1[seurat$anno1 %in% c("EC1", "EC2")] <- 'EC'
unique(seurat$anno1)
unique(seurat$anno1.sub)



seurat_input = seurat
classification = 'anno1'
model_formula = ~  anno1.sub + (1|Patient_ID)
maineffect = 'anno1.sub'
pseudo_factors = c('anno1.sub', 'Patient_ID')



cells.use <- colnames(seurat_input)

object <- seurat_input
labels <- pseudo_factors
assay <- "RNA"
slot <- "counts"

factorlist <- list()
for (i in labels) factorlist[[i]] <- unique(object@meta.data[,i])

meta <- expand.grid(factorlist, stringsAsFactors = FALSE)
rownames(meta) <- apply(meta, 1, function(x) paste0(x, collapse = '.'))

n <- nrow(meta)
out <- matrix(nrow=dim(object[[assay]])[1], ncol=n, data=0)
rownames(out) <- rownames(object[[assay]])
colnames(out) <- rownames(meta)

ncells <- c()
ncounts <- c()
total.cells <- dim(object[[assay]])[2]
for (i in 1:n){
  cells <- 1:total.cells
  for (j in names(meta)) {
    keep  <- which(object@meta.data[[j]] == meta[i,j])
    cells <- cells[cells %in% keep]
  }
  ncells[i] <- length(cells)
  ncounts[i] <- sum(slot(object[[assay]], slot)[,cells])
  
  #some other thing to measure
  if (length(cells)==1) {
    out[,i] <- slot(object[[assay]], slot)[,cells]
  } else {
    out[,i] <- Matrix::rowSums(slot(object[[assay]], slot)[,cells]) # all counts for that gene in this pesudo bulk
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

# don't filter peaks because it is already filtered in archr
keepgenes <- filterByExpr(d$counts, group = d$samples[[maineffect]])
d <- d[keepgenes,]

d <- calcNormFactors(d, method = "TMM")
v <- voomWithDreamWeights(d, model_formula, d$samples, plot=FALSE)
modelfit <- dream(exprObj = v, formula = model_formula, data = d$samples, quiet = TRUE, suppressWarnings = TRUE)
results <- modelfit
rm(d, keepgenes, v, modelfit)



saveRDS(de_results ,'~/gca/yuzhao1/work/final_RC2rna/deg/epithelial/ec/pseudobulk/de_results_POU_EC2vsPOU_EC1.rds')


topTable(de_results[["EC"]],
         coef="anno1.subPOU_EC2",
         p.value = 0.05,
         sort.by = 'logFC',
         number = 10000)




























