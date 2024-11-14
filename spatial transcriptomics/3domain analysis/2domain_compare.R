library(Seurat)
library(IRIS)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')
out.dir <- '~/yuzhao1/work/final_RC2rna/0revision/spatial3/results/'
IRIS_object_list <- readRDS(paste0(out.dir, 'IRIS_object_list.rds'))
seurat_xe <- readRDS('~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged.rds')
seurat_xe$biopsy_location <- ifelse(grepl('pou', seurat_xe$sample), yes = 'pou', no = 'pp')
seurat_xe$cellID <- Cells(seurat_xe)

# specify the group
numCluster <- 5
interested_domain_ID <- '1'
cap_logFC <- 10
cap_padj <- 1e-200

# 1, differential gene for the same domain in different tissues
clusterID <- paste0('nCluster_', numCluster)
IRIS_object <- IRIS_object_list[[clusterID]] 
domain_results <-  IRIS_object@spatialDomain[,c("Slice","spotName","IRIS_domain")]
domain_results$cellID <- paste0(domain_results$Slice, '_', domain_results$spotName)

seurat <- seurat_xe
seurat <- subset(seurat, cells = Cells(seurat)[seurat$cellID %in% domain_results$cellID])
seurat$domain <- mapvalues(seurat$cellID, domain_results$cellID, domain_results$IRIS_domain, warn_missing = F)


seurat <- subset(seurat, cells = Cells(seurat)[seurat$domain == interested_domain_ID])
de_genes <- FindMarkers(object = seurat, ident.1 = 'pou', ident.2 = 'pp', group.by = 'biopsy_location', logfc.threshold = -100000)
de_genes$neg_log10_pval <- -log10(de_genes$p_val)

de_genes$avg_log2FC[de_genes$avg_log2FC > cap_logFC] <- cap_logFC
de_genes$avg_log2FC[de_genes$avg_log2FC < -cap_logFC] <- -cap_logFC
de_genes$p_val_adj[de_genes$p_val_adj < cap_padj] <- cap_padj

# 2. Create the volcano plot
p3 <- EnhancedVolcano(de_genes, rownames(de_genes), labSize = 2.5, 
                      max.overlaps = Inf, drawConnectors=T, arrowheads = F,
                      xlim = c(-cap_logFC, cap_logFC), ylim = c(0, -log10(cap_padj)), 
                      pCutoffCol = 'p_val_adj', pCutoff = 1e-10, 
                      x ="avg_log2FC", 
                      y ="p_val_adj")

pdf(paste0(out.dir, 'domain_deg_volcano/spatial_domain_', numCluster, '_deg_volcano_domain', interested_domain_ID, '.pdf'), width = 8, height = 8)
print(p3)
dev.off()
















