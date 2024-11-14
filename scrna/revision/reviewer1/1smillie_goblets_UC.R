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
seurat <- readRDS('~/spott/yuzhao1/scp259/rds/epi.rds')
seurat <- NormalizeData(seurat)
df_disease_sample_map <- unique(seurat@meta.data[, c('Sample', 'Health')])
rownames(df_disease_sample_map) <- NULL

# build cluster-sample pseudo bulk, filter for pseudo-bulk (now no filter)
seurat$cluster_sample <- paste0(seurat$Cluster, '@', seurat$Sample)
seurat_sub <- subset(seurat, Cluster=='Goblet' & Location=='Epi')
samples_ok <- names(table(seurat_sub$cluster_sample))[table(seurat_sub$cluster_sample) >= 1]
seurat_sub <- subset(seurat_sub, cells=Cells(seurat_sub)[seurat_sub$cluster_sample %in% samples_ok & seurat_sub$Health %in% c("Inflamed","Healthy")])
  
# calculate average expression
labels <- c("ALDOB","REG4","XKR9","TIMP1","SERPINB5","TFF1","AQP3","S100P",
            "LPCAT1","CEACAM6","C2CD4B","TRIM22","B3GALT5","LYZ","TFF2","SCD",
            "TGFBI", "KLK12", "SDR16C5", "PDIA4",  "MYEOV","LINC00261", "B3GNT7","SEC24D","SPINK5",
            "GSTA1","ANPEP","PRAP1","CASD1","RNASE1", "PRDX5","SATB2-AS1","HMGCS2","CA2", "URAD","PPARG")

labels_set1 <- c("ALDOB","REG4","XKR9","TIMP1","SERPINB5","TFF1","AQP3","S100P",
                 "LPCAT1","CEACAM6","C2CD4B","TRIM22","B3GALT5","LYZ","TFF2","SCD",
                 "TGFBI", "KLK12", "SDR16C5", "PDIA4",  "MYEOV","LINC00261", "B3GNT7","SEC24D","SPINK5","CASD1")

labels_set2 <- c("GSTA1","ANPEP","PRAP1", "RNASE1", "PRDX5","SATB2-AS1","HMGCS2","CA2", "URAD","PPARG")

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
Average.expression.mtx <- Average.expression.mtx[labels,]

# create a df for box plots
df <- as.data.frame(Average.expression.mtx)
df$gene <- rownames(df)
melted_df <- melt(df, id='gene')
colnames(melted_df) <- c('gene', 'cluster_sample', 'value')
melted_df$Sample <- melted_df$cluster_sample %>% as.character() %>% strsplit(., split = '@', fixed=T) %>% sapply(.,`[[`,2)
melted_df$Health <- mapvalues(melted_df$Sample, df_disease_sample_map$Sample, df_disease_sample_map$Health, warn_missing = F)

# box plot
color_vec1 <- c("Inflamed"='#F8766D', "Healthy"='#00BFC4')
color_vec2 <- c("Inflamed"='black', "Healthy"='black')
for (i in 1:2) {
  labels_sub <- eval(as.name(paste0('labels_set', i))) %>% sort()
  melted_df_plot <- melted_df[melted_df$gene %in% labels_sub,]
  melted_df_plot$gene <- factor(melted_df_plot$gene, levels = labels_sub)
  pdf(paste0(out.dir, 'Smillie_Goblet_part', i, '.pdf'), height = 6, width = 0.5+length(unique(melted_df_plot$gene))/2)
  p <- ggboxplot(melted_df_plot, x = 'gene', y = 'value', fill = 'Health', outlier.shape = NA, 
            bxp.errorbar=T, bxp.errorbar.width = 0.2)+ 
    scale_fill_manual(values=color_vec1)+
    geom_point(aes(colour = Health), size=0.2,
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





