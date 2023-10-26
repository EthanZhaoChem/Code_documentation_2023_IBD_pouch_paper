dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(Seurat)
library(readr)
library(ggplot2)
library(edgeR)
library(dplyr)
library(plyr)
library(xts)
source('~/yuzhao1/scripts/plot.R')
source('/project/gca/yuzhao1/work/manu/rc2/scripts/umap_colors.R')

df_batch <- read.table('/project/gca/yuzhao1/work/final_RC2rna/metadata/batch.csv', header = T, sep = ' ')
seurat <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/immune3.rds')
seurat$Sample_ID<- paste0(seurat$Patient_ID,'-',seurat$biopsy_location)
seurat$batch <- mapvalues(seurat$Sample_ID, from = df_batch$Sample_ID, to = df_batch$batch)
out.dir <- '/project/gca/yuzhao1/work/final_RC2rna/pca/plots/'

FilePath.Sample_ID.mtx <- paste0(out.dir, 'immune_gene_vs_Sample_ID_RNA_count_mtx.rds') 
FilePath.pdf1 <- paste0(out.dir, 'immune_AllGenes_PC.pdf') 
FilePath.pdf2 <- paste0(out.dir, 'immune_AllGenes_PC1_ByLocation.pdf') 
FilePath.pdf3 <- paste0(out.dir, 'immune_AllGenes_PC2_ByLocation.pdf') 

# egdeR
# creaete empty gene*Sample_ID matrix
nSample_IDs <- length(unique(seurat$Sample_ID))
Sample_ID.mtx <- matrix(ncol = nSample_IDs, nrow = nrow(seurat))
colnames(Sample_ID.mtx) <-  unique(seurat$Sample_ID)
rownames(Sample_ID.mtx) <- rownames(seurat)
ct.mtx <- seurat@assays$RNA@counts
xx <- 0
start.time <- Sys.time()

gene_chunks <- split(rownames(Sample_ID.mtx), ceiling(seq_along(rownames(Sample_ID.mtx))/1000))
for (gene.name in gene_chunks){
  # log
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units = "mins")
  time.taken <- round(time.taken, 2)
  cat(paste0('running time: ', time.taken, ' minutes. ', xx, ' genes out of ', nrow(Sample_ID.mtx), ' has been counted.\n'))
  xx <- xx+length(gene.name)

  # count
  geneQuantity_perCell <- ct.mtx[gene.name,]
  geneQuantity_perCell.t <- t(geneQuantity_perCell)
  geneQuantity_perSample_ID.t <- aggregate(geneQuantity_perCell.t, list(seurat$Sample_ID), sum)
  geneQuantity_perSample_ID <- t(geneQuantity_perSample_ID.t)
  Sample_ID.mtx[gene.name, geneQuantity_perSample_ID['Group.1', ]] <- geneQuantity_perSample_ID[gene.name, ]
}

saveRDS(Sample_ID.mtx, FilePath.Sample_ID.mtx)
Sample_ID.mtx <- readRDS(FilePath.Sample_ID.mtx)

Sample_ID.mtx2 <- matrix(as.numeric(Sample_ID.mtx),ncol = ncol(Sample_ID.mtx))
rownames(Sample_ID.mtx2) <- rownames(Sample_ID.mtx)
colnames(Sample_ID.mtx2) <- colnames(Sample_ID.mtx)

# # if only use variable features
# variable_genes <- VariableFeatures(seurat)
# Sample_ID.mtx2 <- Sample_ID.mtx2[variable_genes, ]

# # if remove MT genes
# MT_genes <- rownames(Sample_ID.mtx2)[grep("^MT-", rownames(Sample_ID.mtx2))]
# Sample_ID.mtx2 <- Sample_ID.mtx2[!(row.names(Sample_ID.mtx2) %in% MT_genes),]

y <- DGEList(counts=Sample_ID.mtx2, group=colnames(Sample_ID.mtx2))
y <- calcNormFactors(y, method = "TMM")
tmm <- cpm(y, normalized.lib.sizes = TRUE, log=T)
res.pca <- prcomp(t(tmm), scale = TRUE)
stdev <- res.pca$sdev
eigValues = (stdev)^2  ## EigenValues
eigValues/sum(eigValues)

# pc1&2 matrix
df <- as.data.frame(res.pca$x)
df$Sample_ID <- rownames(df)
df$patient <- sapply(strsplit(df$Sample_ID, '-'), "[[", 1)
df$biopsy_location <-  sapply(strsplit(df$Sample_ID, '-'), "[[", 2)
df$batch <- mapvalues(df$Sample_ID, df_batch$Sample_ID, df_batch$batch)


pdf(paste0(FilePath.pdf1), width = 5, height = 4)
ggplot(df, aes(x = PC1, y = PC2, label = Sample_ID)) +
  geom_point(aes(shape = batch, fill = biopsy_location), size = 3, alpha = 0.86) +
  scale_fill_manual(values=manual_colors_rc2_location)+
  scale_shape_manual(values = c(21,22,24,25))+
  geom_text_repel(size = 2, max.overlaps = 200)+
  theme_bw()+
  theme(panel.grid =element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 12, face = 'bold'),
        axis.line = element_line(size = 0.5))+
  labs(x = "PC 1 (18.5% explained variance)", y = "PC 2 (10.1% explained variance)", title = "")+
  guides(fill = guide_legend(override.aes = list(shape = c(21))))
dev.off()


df$biopsy_location <- factor(df$biopsy_location, levels = c('AC','POU','PP','TI'))
my_comparisons <- list(c('AC','POU'), c('POU', 'PP'), c('PP', 'TI'))


pdf(paste0(FilePath.pdf2), width = 5, height = 5)
ggboxplot(df, x = "biopsy_location", y = "PC1", fill = "biopsy_location", outlier.shape = NA, 
          bxp.errorbar=T, bxp.errorbar.width = 0.2)+ 
  geom_jitter(color="black", size=1.2, alpha=0.9, width = 0.1) +
  scale_fill_manual(values=manual_colors_rc2_location)+
  stat_compare_means(PC1 ~ biopsy_location, 
                     comparisons = my_comparisons,
                     label.y = c(25, 100, 150),
                     label = "p.signif",                     
                     method = "wilcox.test",
                     paired = F)+
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 10,),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        axis.title = element_text(color = "black", size = 10, angle = 0, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),        legend.position = "none",
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid")
  )+
  labs(x = "Biopsy location", y = "PC 1 (18.5% explained variance)", title = "")
dev.off()

eigValues/sum(eigValues)







