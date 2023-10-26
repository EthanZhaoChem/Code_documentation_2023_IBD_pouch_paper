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
library(VennDiagram)
library(rGREAT)
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
DARs_EC <- readRDS('~/yuzhao1/work/final_RC2atac/peaks/EC_DAR_regions_DifferentContrastLists_ScriptInManuFolder.rds')

################################# analyze peaks per topic #####################
topic_regions_vec <- list()
binary_topic_list <- cisTopicObject@binarized.cisTopics
for (i in seq_along(binary_topic_list)) {
  topic_regions_vec[[i]] <- rownames(binary_topic_list[[i]]) # get the number of rows in the current dataframe and store it in the vector
}

# x2 is a fixed pool to examine x1's falling ratio into it

x2 <- intersect(DARs_EC$EC_POU2vsPP, DARs_EC$EC_ACvsTI) # case colon pouch shared, choose > 20%
x2 <- intersect(DARs_EC$EC_POU1vsAC, DARs_EC$EC_TIvsAC) # case ti and pouch shared, choose > 20%

# colon featrues that are not gained by pouch2: actually didn't plot because these topics are similar
# to the features that gained by pouch (just means they gained but not fully gained)
x2 <- setdiff(DARs_EC$EC_ACvsTI, DARs_EC$EC_POU2vsPP) 
x2 <- setdiff(DARs_EC$EC_TIvsAC, DARs_EC$EC_POU1vsAC) # ileum featrues that are not gained by pouch1
falling_ratio <- list()

for (i in 1:60) {
  a1 <- length(intersect(topic_regions_vec[[i]], x2))
  a2 <- length(topic_regions_vec[[i]])
  falling_ratio[[i]] <- a1/a2
}

# select topics that have a falling ratio >0.2
falling_ratio <- unlist(falling_ratio)
names(falling_ratio) <- paste0('Topic ', 1:60)
falling_ratio %>% sort(., decreasing = T)

# calculate the peaks that fall into the major related topics (as a disection method)
## ileum
x2 <- intersect(DARs_EC$EC_POU1vsAC, DARs_EC$EC_TIvsAC) # case ti and pouch shared, choose > 20%
ileum_topics <- c(45,53,39,60,24,12,46,10,40,20)
DisecctedPeakSet_of_SharedPeaks <- list()
DisecctedPeakSet_of_SharedPeaks_beds <- list()
DisecctedPathways_of_SharedPeaks <- list()
for (i in ileum_topics) {
  gr <- intersect(topic_regions_vec[[i]], x2)
  bed <- data.frame(chr = gr %>%strsplit(., ':') %>% sapply(.,`[[`,1),
                    start = gr %>%strsplit(., ':') %>% sapply(.,`[[`,2) %>%strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric(),
                    end = gr %>%strsplit(., ':') %>% sapply(.,`[[`,2) %>%strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric())
  DisecctedPeakSet_of_SharedPeaks[[as.character(i)]] <- gr
  DisecctedPeakSet_of_SharedPeaks_beds[[as.character(i)]] <- bed
}

saveRDS(DisecctedPeakSet_of_SharedPeaks_beds, '~/yuzhao1/work/final_RC2atac/cistopic/DisecctedPeakSet_of_SharedPeaks_beds.rds')

DisecctedPeakSet_of_SharedPeaks_beds <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/DisecctedPeakSet_of_SharedPeaks_beds.rds')
for (i in ileum_topics) {
  bed <- DisecctedPeakSet_of_SharedPeaks_beds[[as.character(i)]]
  job = submitGreatJob(bed, species = "hg38", rule = "oneClosest")
  tb = getEnrichmentTables(job)
  pathways <- tb[["GO Biological Process"]]
  pathways <- pathways[pathways$Hyper_Adjp_BH<0.01,]
  pathways <- pathways[pathways$Hyper_Fold_Enrichment>1,]
  DisecctedPathways_of_SharedPeaks[[as.character(i)]] <- pathways
}
saveRDS(DisecctedPathways_of_SharedPeaks, '~/yuzhao1/work/final_RC2atac/cistopic/DisecctedPathways_of_SharedPeaks_ileum.rds')


explained_ratio <- length(unique(unlist(DisecctedPeakSet_of_SharedPeaks)))/length(x2)






################################# shared topics rGREAT results #########################
pdf(paste0(plot.dir, 'great_colon.pdf'), height = 8.5, width = 8)
ontologyDotPlot(cisTopicObject, top=5, topics=c(23,42,30,28,1,14,19), var.y='name', order.by='Binom_Adjp_BH',
                min.size = 2,
                max.size = 10)
dev.off()

pdf(paste0(plot.dir, 'great_heatmapBasedShared_colon.pdf'), height = 8.5, width = 15)
ontologyDotPlot(cisTopicObject, top=5, topics=c(7, 2, 30, 52, 27, 43, 34, 48, 59, 14, 31, 28, 6, 55), var.y='name', order.by='Binom_Adjp_BH',
                min.size = 2,
                max.size = 10)
dev.off()

pdf(paste0(plot.dir, 'great_heatmapBasedShared_ileum.pdf'), height = 8.5, width = 10)
ontologyDotPlot(cisTopicObject, top=5, topics=c(25, 40, 45, 46, 53, 10, 3), var.y='name', order.by='Binom_Adjp_BH',
                min.size = 2,
                max.size = 10)
dev.off()

pdf(paste0(plot.dir, 'great_ileum.pdf'), height = 8.5, width = 12)
ontologyDotPlot(cisTopicObject, top=5, topics=c(45,53,39,60,24,12,46,10,40,20), var.y='name', order.by='Binom_Adjp_BH',
                min.size = 2,
                max.size = 6)
dev.off()

################################# DAR enrichment rGREAT results #########################
# these peaks are too many to generate a specific GO process
gr <- DARs_EC$EC_ACvsTI # case colon 
gr <- DARs_EC$EC_TIvsAC # case ileum 
gr <- intersect(DARs_EC$EC_POU2vsPP, DARs_EC$EC_ACvsTI) # case colon pouch shared
gr <- intersect(DARs_EC$EC_POU1vsAC, DARs_EC$EC_TIvsAC) # ileum pouch shared

bed <- data.frame(chr = gr %>%strsplit(., ':') %>% sapply(.,`[[`,1),
                 start = gr %>%strsplit(., ':') %>% sapply(.,`[[`,2) %>%strsplit(., '-') %>% sapply(.,`[[`,1) %>% as.numeric(),
                 end = gr %>%strsplit(., ':') %>% sapply(.,`[[`,2) %>%strsplit(., '-') %>% sapply(.,`[[`,2) %>% as.numeric())

job = submitGreatJob(bed, species = "hg38", rule = "oneClosest")
tb = getEnrichmentTables(job)

tb[["GO Biological Process"]] %>% View()
res = plotRegionGeneAssociationGraphs(job)

saveRDS(tb, '~/yuzhao1/work/final_RC2atac/cistopic/enrichment_ileum_pouch.rds')

tb <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/enrichment_colon_pouch.rds')
tb <- readRDS('~/yuzhao1/work/final_RC2atac/cistopic/enrichment_ileum_pouch.rds')
xx <- tb[["GO Biological Process"]]

xx <- xx[xx$Hyper_Adjp_BH<0.01,]
xx <- xx[xx$Hyper_Fold_Enrichment>2,]
View(xx)
























