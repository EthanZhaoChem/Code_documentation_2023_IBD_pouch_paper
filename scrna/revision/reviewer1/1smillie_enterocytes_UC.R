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
seurat_ec <- subset(seurat, Cluster=='Enterocytes' & Location=='Epi')
samples_ok <- names(table(seurat_ec$cluster_sample))[table(seurat_ec$cluster_sample) >= 1]
seurat_ec <- subset(seurat_ec, cells=Cells(seurat_ec)[seurat_ec$cluster_sample %in% samples_ok & seurat_ec$Health %in% c("Inflamed","Healthy")])
  
# calculate average expression
labels <- c(
  'REG4','DMBT1','SLC6A14','VNN1','TGM2','TFF1','HLA-DRB1','MGAT3','PI3','OLFM4','LPCAT1',
  'UNC5CL','C2','TRIM22','C4BPB','LCN2','MYEOV','CASP1','UNC13D','CFB','TRIM40','SESTD1','KDELR3',
  'LPIN1','ARFGAP3','TMEM92','S100A11', 'ELOVL7','GBP1','FFAR4','ITGA2','XDH','DUSP6','EFNA2','STAT1',
  'MAP3K5','SLC20A1','ZNF703','P2RX4','RNF186','MESP1','MUC4','SELENBP1','OASL','DEFB1','HMGCS2',
  'CA2','AQP8','ANPEP','CASD1','SPINK5','KAZN','SMAD9','HYI','GPX4','TNNC2')

labels_set1 <- c('SMAD9','KAZN','SPINK5','MUC4','CASD1','MAP3K5','STAT1','EFNA2',
                 'ITGA2','FFAR4','LPIN1','UNC13D','CASP1','MYEOV','OLFM4',
                 'PI3','TFF1','SLC6A14','DUSP6','S100A11')

labels_set2 <- c('LCN2','C4BPB','XDH','GBP1','ELOVL7','TMEM92','ARFGAP3',
                 'KDELR3','SESTD1','TRIM40','CFB','TRIM22','C2','UNC5CL',
                 'LPCAT1','MGAT3','HLA−DRB1','TGM2','VNN1','DMBT1','REG4')

labels_set3 <- c('P2RX4','SLC20A1','AQP8','CA2','HMGCS2',
                 'DEFB1','OASL','SELENBP1','MESP1','RNF186',
                 'ZNF703','TNNC2','GPX4','HYI','ANPEP')

Average.expression.mtx <- AverageExpression(
  seurat_ec,
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
for (i in 1:3) {
  labels_sub <- eval(as.name(paste0('labels_set', i))) %>% sort()
  melted_df_plot <- melted_df[melted_df$gene %in% labels_sub,]
  melted_df_plot$gene <- factor(melted_df_plot$gene, levels = labels_sub)
  pdf(paste0(out.dir, 'Smillie_EC_part', i, '.pdf'), height = 6, width = 0.5+length(unique(melted_df_plot$gene))/2)
  p <- ggboxplot(melted_df_plot, x = 'gene', y = 'value', fill = 'Health', outlier.shape = NA, 
            bxp.errorbar=T, bxp.errorbar.width = 0.2)+ 
    scale_fill_manual(values=color_vec1)+
    geom_point(aes(colour = Health), size=0.2,
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2, jitter.height = 0.1,))+
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





