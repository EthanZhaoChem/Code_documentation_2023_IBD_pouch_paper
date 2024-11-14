library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

source('~/yuzhao1/work/final_RC2rna/0revision/spatial2/0helper.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')
source('~/yuzhao1/scripts/plot.R')

plot.dir <- '~/yuzhao1/work/final_RC2rna/0revision/spatial2/plots/'
seurat <- readRDS('~/yuzhao1/work/final_RC2rna/0revision/spatial2/rds/xenium.merged_processed.rds')

# 1. group EC
seurat$cluster <- seurat$predicted.celltype
seurat$cluster[seurat$cluster %in% c('EC1-1','EC1-2', 'EC2-1', 'EC2-2')] <- 'EC'
seurat.ec <- subset(seurat, cluster=='EC')
seurat.ec@meta.data[seurat.ec$predicted.celltype %in% c('EC1-1','EC1-2'), 'cluster'] <- 'POU-EC1'
seurat.ec@meta.data[seurat.ec$predicted.celltype %in% c('EC2-1','EC2-2'), 'cluster'] <- 'POU-EC2'
seurat.ec@meta.data[seurat.ec$sample == 'ac', 'cluster'] <- 'AC-EC'
seurat.ec@meta.data[seurat.ec$sample == 'ti', 'cluster'] <- 'TI-EC'
seurat.ec@meta.data[seurat.ec$sample %in% c("ppR1", "ppR2"), 'cluster'] <- 'PP-EC'

Pgenes <- c('DMBT1',  'MUC2', 'OLFM4',  'LCN2',
            'ADH1C', 'SATB2', 'CEACAM7', 'CA2', 'AQP8', 'HES1', 'NXPE4', 'HMGCS2',
            'RPSA', 'CFTR', 'RPL11', 'CEACAM5', 'SOX9', 'MECOM', 'TFCP2L1', 'CD24','FCGBP',
            'GSTA1', 'ACE2', 'RBP2', 'APOA4', 'FABP6', 'HLA−DRA','ABCC2', 'ENPEP',
            'SELENOP', 'TRPM6', 
            "REG1B", "FOXD1", "GATA5", "KLF7", "REG4", "IRF8", 
            "STAT6", "STAT5", "SATB1", 
            "HOXB7","MT1A", "KLF11" ,"MT1A", "ABCB9",
            "REG3A","TUBA1B",  "CELA3B",
            "FABP1", "KRT20","ZBTB7B",
            "DEFA6", "DEFA5", "CLCA1", "TFF1", 'ATOH1',
            'GUCA2B') %>% intersect(rownames(seurat))

# 2. plot vln plot
p_vln <- list()
for (Pgenes.i in Pgenes) {
  p_vln[[Pgenes.i]] <- VlnPlot(seurat.ec, Pgenes.i, group.by = 'sample')
}
pdf(paste0(plot.dir, 'vln', '.pdf'), width = 6, height = 5)
print(p_vln)
dev.off()


# 3. plot heatmap

Average.expression.mtx <- AverageExpression(
  seurat.ec,
  assays = 'RNA',
  features = Pgenes,
  return.seurat = FALSE,
  group.by = "cluster",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)
Average.expression.mtx <- Average.expression.mtx$RNA %>% t(.) %>% scale(.) %>% t(.)
labels.idx <- match(Pgenes, rownames(Average.expression.mtx))

pdf(paste0(plot.dir, 'heatmap.pdf'), width = 5, height = 8)
p1 <- Heatmap(Average.expression.mtx,
              col = rc2_rna_heatmap_colors_gradient1,
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, show_row_names = T,
              row_title = NULL,
              row_gap = unit(1, "mm"), column_gap = unit(0, "mm"), 
              border = T,
              border_gp = gpar(col = "black", lty = 1, lwd = 0),
              
              heatmap_legend_param = list(title = "Mean expression z score"),
              use_raster = F)
print(p1)
dev.off()

