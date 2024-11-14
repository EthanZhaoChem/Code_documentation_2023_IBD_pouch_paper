dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(harmony)
library(Seurat)
library(SeuratObject)
library(reshape2)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/helper_seurat.R')
source('~/yuzhao1/scripts/seurat/deg_pseudobulk.R')
source('/project/gca/yuzhao1/work/manu/rc2/scripts/rna_annotation_markers.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')

out.dir <- '~/yuzhao1/work/final_RC2rna/0revision/reviewer1/plots/'

# read and normalize data
seurat <- readRDS('~/spott/yuzhao1/broadSCP1884/dataset/epithelial.rds')
seurat <- NormalizeData(seurat)
seurat$Cluster <- seurat$anno1
seurat$Sample <- seurat$Sample_ID
seurat$cluster_sample <- paste0(seurat$Cluster, '@', seurat$Sample)
df_location_sample_map <- unique(seurat@meta.data[, c('biopsy_location', 'Sample')])
rownames(df_location_sample_map) <- NULL

# build cluster-sample pseudo bulk, filter for pseudo-bulk (now no filter)
ec_clusters <- c("Enterocytes CA1 CA2 CA4-", "Enterocytes TMIGD1 MEP1A", "Enterocytes TMIGD1 MEP1A GSTA1")
seurat_sub <- subset(seurat, cells=Cells(seurat)[seurat$anno1 %in% ec_clusters & seurat$disease_status %in% c("normal") & seurat$Layer=='E'])
samples_ok <- names(table(seurat_sub$cluster_sample))[table(seurat_sub$cluster_sample) >= 1]
seurat_sub <- subset(seurat_sub, cells=Cells(seurat_sub)[seurat_sub$cluster_sample %in% samples_ok])
  
# calculate average expression
labels <- c('STAT6','IRF8','KLF7','GATA5','FOXD1','REG1B','FABP6','GUCA2B',
            'TRPM6','SELENOP','ACE2','GSTA1','ENPEP','ABCC2','RBP2','APOA4',
            'FABP1','SATB2','HES1','LCN2','CA2','NXPE4','CEACAM7','AQP8','SOX9','RPSA','CFTR',
            'MECOM','ADH1C','FCGBP','HMGCS2','CD24','CEACAM5','TFCP2L1')

labels_set1 <- c('STAT6','IRF8','KLF7','GATA5','FOXD1','REG1B','FABP6','GUCA2B',
                 'TRPM6','SELENOP','ACE2','GSTA1','ENPEP','ABCC2','RBP2','APOA4')

labels_set2 <- c('FABP1','SATB2','HES1','LCN2','CA2','NXPE4','CEACAM7','AQP8','SOX9','RPSA','CFTR',
                 'MECOM','ADH1C','FCGBP','HMGCS2','CD24','CEACAM5','TFCP2L1')

Average.expression.mtx <- AverageExpression(
  seurat_sub,
  assays = 'RNA',
  return.seurat = FALSE,
  group.by = "cluster_sample",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)
Average.expression.mtx <- Average.expression.mtx$RNA %>% t(.) %>% scale(.) %>% t(.)
labels <- intersect(rownames(Average.expression.mtx), labels)
labels_set1 <- intersect(rownames(Average.expression.mtx), labels_set1)
labels_set2 <- intersect(rownames(Average.expression.mtx), labels_set2)
Average.expression.mtx <- Average.expression.mtx[labels,]

# create a df for box plots
df <- as.data.frame(Average.expression.mtx)
df$gene <- rownames(df)
melted_df <- melt(df, id='gene')
colnames(melted_df) <- c('gene', 'cluster_sample', 'value')
melted_df$Sample <- melted_df$cluster_sample %>% as.character() %>% strsplit(., split = '@', fixed=T) %>% sapply(.,`[[`,2)
melted_df$biopsy_location <- mapvalues(melted_df$Sample, df_location_sample_map$Sample, df_location_sample_map$biopsy_location, warn_missing = F)

# box plot
color_vec1 <- c("colon"='#F8766D', "ileum"='#00BFC4')
color_vec2 <- c("colon"='black', "ileum"='black')
for (i in 1:2) {
  labels_sub <- eval(as.name(paste0('labels_set', i))) %>% sort()
  melted_df_plot <- melted_df[melted_df$gene %in% labels_sub,]
  melted_df_plot$gene <- factor(melted_df_plot$gene, levels = labels_sub)
  pdf(paste0(out.dir, 'Kong_EC_part', i, '.pdf'), height = 6, width = 0.5+length(unique(melted_df_plot$gene))/2)
  p <- ggboxplot(melted_df_plot, x = 'gene', y = 'value', fill = 'biopsy_location', outlier.shape = NA, 
            bxp.errorbar=T, bxp.errorbar.width = 0.2)+ 
    scale_fill_manual(values=color_vec1)+
    geom_point(aes(colour = biopsy_location), size=0.2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0.1))+
    scale_color_manual(values=color_vec2)+
    theme_pubr()+
    theme(axis.text.y = element_text(size=8),
          axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
          axis.title.y = element_text(size=8),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "bottom",
          plot.title = element_text(size=20, hjust=0.5, face = 'bold'))+
    labs(title = paste0(), y='Normalized expression (Z score)')
  print(p)
  dev.off()
}





