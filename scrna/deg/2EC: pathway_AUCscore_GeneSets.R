dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')

packages <- c('dplyr', 'Seurat', 'patchwork', 'ggplot2','ggpubr', 'ggsci', 'plyr', 'stringr' )
library(AUCell)
library(GSEABase)
library(pheatmap)
library(ggrepel)
library(Matrix)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
lapply(packages, library, character.only = TRUE)
source('~/yuzhao1/scripts/plot.R')

# 0. helper functions
ethan_get_auc_score <- function(cells_ranking, geneset, aucMaxRank_manual, geneset_name = 'TBD'){
  geneset <- unique(geneset)
  geneSet <- GeneSet(geneset, setName = geneset_name)
  cells_AUC <- AUCell_calcAUC(geneSet, cells_ranking, aucMaxRank = ceiling(aucMaxRank_manual * nrow(cells_ranking)))
  AUC_scores <- getAUC(cells_AUC) 
  score_temp <- as.data.frame(t(AUC_scores))
  colnames(score_temp) <- c(geneset_name)
  results <- list()
  results[[1]] <- score_temp
  results[[2]] <- cells_AUC
  return(results)
}

# 1. save/read msigdb.human to resource folder
# msigdb.human <- msigdbr(species = "Homo sapiens")
# saveRDS(msigdb.human, '/project/gca/yuzhao1/resource/msigdb/msigdb.human.rds')
union <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
msigdb.human <- readRDS('/project/gca/yuzhao1/resource/msigdb/msigdb.human.rds')
cells_ranking <-readRDS('~/yuzhao1/work/final_RC2rna/deg/union/union_cells_ranking.rds')

seurat <- union
# 2. extract gene names by pathway names
genesets <- list()
pathway_name <- 'HALLMARK_INFLAMMATORY_RESPONSE'
genesets[[pathway_name]] <- msigdb.human [msigdb.human $gs_name %in% c(pathway_name), ]$gene_symbol

# 3. score a pathway
pathway_name <- 'HALLMARK_INFLAMMATORY_RESPONSE'
results <- ethan_get_auc_score(cells_ranking = cells_ranking,
                               geneset = genesets[[pathway_name]],
                               aucMaxRank_manual = 0.05,
                               geneset_name  = pathway_name)

scoreDf <- results[[1]]

seurat <- AddMetaData(
  object = seurat,
  metadata = scoreDf,
  col.name = names(scoreDf)
)
















