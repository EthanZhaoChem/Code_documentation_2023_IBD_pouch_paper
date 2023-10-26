dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')

packages <- c('dplyr', 'Seurat', 'patchwork', 'ggplot2','ggpubr', 'ggsci', 'plyr', 'stringr' )
library(AUCell)
library(GSEABase)
library(pheatmap)
library(ggrepel)
library(Matrix)
lapply(packages, library, character.only = TRUE)
source('~/yuzhao1/scripts/plot.R')

union <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
out.dir <- '~/yuzhao1/work/final_RC2rna/deg/union/'
seurat <- union

norm_counts <- seurat@assays$RNA@counts
xx <- as.matrix(norm_counts[, c(1:10000)])
NN <- ceiling(ncol(norm_counts)/10000)
for (ii in 2:NN){
  nstart <- (ii-1)*10000+1
  nstop <- min(ncol(norm_counts), ii*10000)
  xx <- cbind(xx,as.matrix(norm_counts[, nstart:nstop]))
}
exprMatrix <- xx
rm(xx, ii, NN, nstart, nstop)

# for ranking, rounding up is accecptable
exprMatrix.rounded <- round(exprMatrix, digits = 0)
cells_ranking <- AUCell_buildRankings(exprMatrix.rounded, plotStats = T)



saveRDS(cells_ranking, paste0(out.dir, 'union_cells_ranking.rds'))
# cells_ranking <-readRDS('union_cells_ranking.rds')


