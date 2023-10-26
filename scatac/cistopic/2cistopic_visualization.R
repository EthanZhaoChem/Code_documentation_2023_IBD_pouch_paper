library(Matrix)
library(cisTopic)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(plyr)
library(stringr)
library(ArchR)
library(R.utils)
library(rtracklayer)
library(AUCell)
library(RcisTarget)
source('~/yuzhao1/scripts/plot.R')

# prepare
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
plot.dir <- '/project/gca/yuzhao1/work/final_RC2atac/cistopic/plots/'
cisTopicObject <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic60_enriched_great_cistrome.rds') # default model60
# cisTopicObject@selected.model <- cisTopicObject@models$`60`
dic_cellname_to_anno1 <- list(
  cellnames = rownames(proj@cellColData),
  anno1s = as.data.frame(proj@cellColData)$anno1
)

################################# heatmap plot ################################
# one to one linked cellname and anno1
cellnames_subset <- proj$cellNames
anno1s <- dic_cellname_to_anno1$anno1s[match(cellnames_subset, dic_cellname_to_anno1$cellnames)]
metadata_ordered <- as.data.frame(proj@cellColData)[match(cellnames_subset, dic_cellname_to_anno1$cellnames),]

# make cell topic score mtx
idx_cis <- match(cellnames_subset, colnames(cisTopicObject@selected.model$document_sums))
mtx <- cisTopicObject@selected.model$document_sums[,idx_cis]
rownames(mtx) <- paste0('Topic_', 1:nrow(cisTopicObject@selected.model$topics))

# use anno1 as the 'variable name' for the group you wanna use
mtx_cellrow <- as.data.frame(t(mtx))
mtx_cellrow$anno1 <- metadata_ordered$anno1.loc
mtx_cellrow <- mtx_cellrow[which(mtx_cellrow$anno1 %in% c("EC-POU1", "EC-POU2", "EC-TI", "EC-PP", "EC-AC")), ]

mtx_anno1row <- mtx_cellrow %>%
  group_by(anno1) %>%
  summarise_all(funs(sum))
mtx_anno1row <- as.data.frame(mtx_anno1row)
rownames(mtx_anno1row) <- mtx_anno1row$anno1
mtx_anno1row$anno1 <- NULL

# regress out the clusuter original size == use proportion
for (celltype.i in rownames(mtx_anno1row)){
  cluster.size <- sum(mtx_cellrow$anno1 == celltype.i)
  mtx_anno1row[celltype.i,] <- as.numeric(mtx_anno1row[celltype.i,])/cluster.size
}
mtx_scaled <- t(scale(mtx_anno1row)) # column-wise scale, scale in topic
mtx_scaled2 <- scale(t(mtx_anno1row)) # scale in celltype

pdf(paste0(plot.dir, 'EC_60topic_scaled_in_topic.pdf'), height = 12, width = 4)
col_fun = colorRamp2(-8:8/3, manual_colors_gradient1)
p <- Heatmap(mtx_scaled, cluster_rows = T, show_row_dend = T, cluster_columns = T, name = 'topic', col = col_fun)
print(p)
dev.off()

pdf(paste0(plot.dir, 'EC_60topic_scaled_in_celltype.pdf'), height = 12, width = 4)
col_fun = colorRamp2(-8:8/3, manual_colors_gradient1)
p <- Heatmap(mtx_scaled2, cluster_rows = T, show_row_dend = T, cluster_columns = T, name = 'topic', col = col_fun)
print(p)
dev.off()


################################# read TF results #############################


DT::datatable(
  cisTopicObject@binarized.RcisTarget[[45]][, -c("enrichedRegions", "TF_lowConf"), with =
                                              FALSE],
  escape = FALSE,
  filter = "top",
  options = list(pageLength = 100)
) 


################################# read cistrome results #######################
aucellRankings <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/aucellRankings.rds')


cisTopicObject <- getCistromeEnrichment(cisTopicObject, topic=30, TFname='JUN', aucellRankings = aucellRankings, aucMaxRank = 0.05*nrow(aucellRankings), plot=FALSE)

xx <- colnames(cisTopicObject@cell.data)[ncol(cisTopicObject@cell.data)]

pdf(paste0(plot.dir, xx, '_cistrome.pdf'), height = 8.5, width = 12)
df <- data.frame((proj@cellColData))
df$biopsy_location <- factor(proj$biopsy_location, levels = c("AC", "TI", "POU", "PP"))
df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
# retrive data from the correct nTopic
idx <- match(proj$cellNames, rownames(cisTopicObject@cell.data))
df$feature_to_plot <- cisTopicObject@cell.data[idx, xx]
p <- plot_df_umap_custom(df, plot_feature = T)+
  facet_wrap( ~ biopsy_location , ncol=2) +
  theme(strip.background = element_rect(fill = "white", colour = "white"), strip.text = element_text(size = 12))
print(p)
dev.off()


################################# analyze peaks per topic #####################
xx <- cisTopicObject@binarized.cisTopics$Topic59


















