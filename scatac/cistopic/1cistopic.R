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
source('~/yuzhao1/scripts/helper_archr.R')

# has been run based on 2022-10-06 peakset

proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
plot.dir <- '/project/gca/yuzhao1/work/final_RC2atac/cistopic/plots/'

# # # ################################## get unfiltered_peak.mtx ##############################################
# addArchRThreads(1)
# peak.mtx <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = T)
# saveRDS(peak.mtx, paste0(proj@projectMetadata$outputDirectory, '/peakMat_Anno2_binarized.rds'))


# # ################################## get unfiltered_peak.mtx ##############################################
# # cistopic requires binarized matrix
# peak.mtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peakMat_Anno2_binarized.rds'))
# counts <- peak.mtx@assays@data$PeakMatrix
# rownames(counts) <- as.character(peak.mtx@rowRanges, ignore.strand=FALSE)
# cisTopicObject <- createcisTopicObject(counts, project.name='gut')
# gca_metadata <- as.data.frame(proj@cellColData)
# cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = gca_metadata)
# cisTopicObject <- runWarpLDAModels(cisTopicObject,
#                                    topic=c(12, 15, 18),
#                                    seed=6,
#                                    nCores=1,
#                                    iterations = 500,
#                                    addModels=FALSE,
#                                    tmp = '~/yuzhao1/work/final_RC2atac/cistopic/models/')
# saveRDS(cisTopicObject, '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_12_15_18.rds')



################################## add metadata ##############################################
# cisTopicObject <- readRDS( '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject.rds')
# cisTopicObject@selected.model <- cisTopicObject@models$`20`
# addinfo <- data.frame(row.names = proj$cellNames, predictedCell_Co_sample = proj$predictedGroup_Co_sample)
# cisTopicObject <- addCellMetadata(cisTopicObject, addinfo)

# ################################## assess and select models ##############################################
# # select 60 after comparing model perplexity
# cisTopicObject <- readRDS( '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject.rds')
# for (xx in c(10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500)){
#   xx <- as.character(xx)
#   cat("\n")
#   cat(xx, '\n')
#   cat(cisTopicObject@models[[xx]]$perplexity)
#   cat("\n")
# }
# 
# cisTopicObject@selected.model <- cisTopicObject@models$`18`
# plots <- list()
# for (topic_id in 1:nrow(cisTopicObject@selected.model$topics)){
#   df <- data.frame((proj@cellColData))
#   df$biopsy_location <- factor(proj$biopsy_location, levels = c("AC", "TI", "POU", "PP"))
#   df$embedding1 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_1`
#   df$embedding2 <- proj@embeddings$Harmony_UMAP$df$`Harmony#UMAP_Dimension_2`
#   # retrive data from the correct nTopic
#   idx <- match(proj$cellNames, colnames(cisTopicObject@selected.model$document_sums))
#   df$feature_to_plot <- cisTopicObject@selected.model$document_sums[topic_id, idx]
#   p <- plot_df_umap_custom(df, plot_feature = T)+
#     facet_wrap( ~ biopsy_location , ncol=2) +
#     theme(strip.background = element_rect(fill = "white", colour = "white"), strip.text = element_text(size = 12))
#   plots[[topic_id]] <- p
# }
# pdf(paste0('/project/gca/yuzhao1/work/final_RC2atac/cistopic/topic18/', 'TopicUmapByLoc', '.pdf'), height = 8.5, width = 7)
# print(plots)
# dev.off()

# ################################## structure plot ##############################################
# cisTopicObject@selected.model <- cisTopicObject@models$`20`
# idx <- which(cisTopicObject@cell.data[['anno1_loc']] %in% c("EC-POU1", "EC-POU2", "EC-TI", "EC-PP", "EC-AC"))
# df <- as.data.frame(t(cisTopicObject@selected.model$document_sums[, idx]))
# colnames(df) <- paste0('topic', 1:20)
# xx <- cisTopicObject@cell.data[['anno1_loc']][idx]
# df$anno1_loc <- xx
# saveRDS(df, '~/yuzhao1/work/final_RC2atac/cistopic/structure/df_topic20.rds')

# # ################################## calculate umap based on topics ##############################################
# cisTopicObject <- readRDS( '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject.rds')
# cisTopicObject@selected.model <- cisTopicObject@models$`60`
# cisTopicObject <- runUmap(cisTopicObject, target='cell')
# 
# png(paste0(plot.dir,'umap_location.png'), res = 300, height = 3000, width = 3000)
# df <- data.frame((cisTopicObject@cell.data))
# df$biopsy_location <- factor(df$biopsy_location, levels = c("AC", "TI", "POU", "PP"))
# df$embedding1 <- cisTopicObject@dr$cell$Umap[,1]
# df$embedding2 <- cisTopicObject@dr$cell$Umap[,2]
# df$cluster_name <- df$anno1
# p <- plot_df_umap_custom(df, show.label = 'name')+
#   facet_wrap( ~ biopsy_location , ncol=2) +
#   theme(strip.background = element_rect(fill = "white", colour = "white"), strip.text = element_text(size = 12))
# print(p)
# dev.off()

# ################################## plot cell topic heat map ##############################################
# cisTopicObject <- readRDS( '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject.rds')
# cisTopicObject@selected.model <- cisTopicObject@models$`60`
# 
# dic_cellname_to_anno1 <- list(
#   cellnames = rownames(proj@cellColData),
#   anno1s = as.data.frame(proj@cellColData)$anno1
# )
# cellnames_subset <- proj$cellNames
# anno1s <- dic_cellname_to_anno1$anno1s[match(cellnames_subset, dic_cellname_to_anno1$cellnames)]
# idx_cis <- match(cellnames_subset, colnames(cisTopicObject@selected.model$document_sums))
# mtx <- cisTopicObject@selected.model$document_sums[,idx_cis]
# rownames(mtx) <- paste0('Topic_', 1:nrow(cisTopicObject@selected.model$topics))
# 
# mtx_cellrow <- as.data.frame(t(mtx))
# mtx_cellrow$anno1 <- anno1s
# mtx_anno1row <- mtx_cellrow %>%
#   group_by(anno1) %>%
#   summarise_all(funs(sum))
# mtx_anno1row <- as.data.frame(mtx_anno1row)
# rownames(mtx_anno1row) <- mtx_anno1row$anno1
# mtx_anno1row$anno1 <- NULL
# 
# # regress out the clusuter original size == use proportion
# for (celltype.i in rownames(mtx_anno1row)){
#   cluster.size <- sum(proj$anno1 == celltype.i)
#   mtx_anno1row[celltype.i,] <- as.numeric(mtx_anno1row[celltype.i,])/cluster.size
# }
# mtx_scaled <- t(scale(mtx_anno1row)) # column-wise scale, scale in topic
# mtx_scaled2 <- t(scale(t(mtx_anno1row))) # scale in celltype
# 
# png(paste0(plot.dir, '60topic_scaled_in_topic.png'), res = 300, height = 3600, width = 3000)
# col_fun = colorRamp2(-8:8/3, manual_colors_gradient1)
# p <- Heatmap(mtx_scaled, cluster_rows = T, show_row_dend = T, cluster_columns = T, name = 'topic', col = col_fun)
# print(p)
# dev.off()
# 
# png(paste0(plot.dir, '60topic_scaled_in_celltype.png'), res = 300, height = 2600, width = 6000)
# col_fun = colorRamp2(-8:8/3, manual_colors_gradient1)
# p <- Heatmap(mtx_scaled2, cluster_rows = T, show_row_dend = T, cluster_columns = T, name = 'topic', col = col_fun)
# print(p)
# dev.off()


################################## Motif analysis ##############################################
cisTopicObject <- readRDS( '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject.rds')
cisTopicObject@selected.model <- cisTopicObject@models$`20`

# basic steps for Analysis of the regulatory topics
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=F)

# Import chain file
hg38ToHg19.chain <- import.chain('~/yuzhao1/resource/hg38ToHg19.over.chain')


# Obtain liftOver dictionary (as list)
hg38_coord <- cisTopicObject@region.ranges
hg38_to_hg19_list <- liftOver(hg38_coord, hg38ToHg19.chain)

# Run cisTopic based on liftover to hg19 coordinates
cisTopicObject <- binarizedcisTopicsToCtx(cisTopicObject, liftOver=hg38_to_hg19_list, genome="hg19")
cisTopicObject <- scoredRegionsToCtx(cisTopicObject, liftOver=hg38_to_hg19_list, genome="hg19")
pathToFeather <- "~/yuzhao1/resource/hg19-regions-9species.all_regions.mc9nr.feather"

# cisTopicObject,code chunk in core=1 is deprecated, use ncores>1
cisTopicObject <- topicsRcisTarget(cisTopicObject, genome='hg19', pathToFeather,
                                   reduced_database=FALSE, nesThreshold=3, rocthr=0.005, maxRank=20000, nCores=1)
cisTopicObject<- getCistromes(cisTopicObject, annotation = 'Both', nCores=1)
saveRDS(cisTopicObject,
        '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic20_enriched.rds')



################################# read TF results ##############################################
#
cisTopicObject <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic20_enriched.rds')

DT::datatable(
  cisTopicObject@binarized.RcisTarget[[15]][, -c("enrichedRegions", "TF_lowConf"), with =
                                              FALSE],
  escape = FALSE,
  filter = "top",
  options = list(pageLength = 100)
)

################################# TF enrichment in archr ##############################################
cisTopicObject <- readRDS( '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic20_enriched.rds')
topic_enrichment_results <- list()
for (single_topic in paste0('Topic', 1:20)){
  topic_region_vector <- rownames(cisTopicObject@binarized.cisTopics[[single_topic]])
  topic_enrichment <- archr_customized_motif_enrichment(ArchRProj = proj,
                                                        peakAnnotation = "Motif",
                                                        candidate_peaks_vector = topic_region_vector)
  topic_enrichment_results[[single_topic]] <- topic_enrichment
}
saveRDS(topic_enrichment_results, '~/yuzhao1/work/final_RC2atac/cistopic/topic20/TF_topic_enrichment_results.rds')

################################# read rGREAT results ##############################################
# # # Run GREAT based on liftover to hg19 coordinates
# cisTopicObject <- GREAT(cisTopicObject, genome='hg19', liftOver=hg38_to_hg19_list, fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
# saveRDS(cisTopicObject,
#         '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic60_enriched_great.rds')

# # read great results
# cisTopicObject <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic60_enriched_great.rds')
# 
# 
# pdf(paste0(plot.dir, 'great.pdf'), height = 8.5, width = 12)
# ontologyDotPlot(cisTopicObject, top=5, topics=c(30, 36, 59), var.y='name', order.by='Binom_Adjp_BH',
#                 min.size = 2,
#                 max.size = 6)
# dev.off()

# ################################# run cistrome ##############################################
# library(AUCell)
# cisTopicObject <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic60_enriched_great.rds')
# # Cistrome defines a cistrome as a set of sequences enriched for motifs linked to a certain transcription factor.
# cisTopicObject <- getCistromes(cisTopicObject, annotation = 'Both', nCores=5)
# 
# # Compute AUC rankings based on the predictive distribution
# pred.matrix <- predictiveDistribution(cisTopicObject)
# aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE, nCores = 1)
# 
# # save aucellRankings and cisTopicObject
# saveRDS(cisTopicObject,
#         '~/yuzhao1/work/final_RC2atac/cistopic/cisTopicObject_topic60_enriched_great_cistrome.rds')
# saveRDS(aucellRankings, '~/yuzhao1/work/final_RC2atac/cistopic/aucellRankings.rds')



























