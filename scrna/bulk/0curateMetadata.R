dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')

################ curate old metadata ##############
# metadata <- read.table('/project/gca/yuzhao1/work/final_RC2rna/bulk/metadata/GSE81266_series_matrix.txt',
#                        sep = '\t', header = F, comment.char = '#') %>% t(.)
# colnames(metadata) <- metadata[1,]
# colnames(metadata) <- gsub('!', '', colnames(metadata)) 
# colnames(metadata) <- gsub('_ch1', '', colnames(metadata)) 
# metadata <- metadata[-1,]
# colnames(metadata)[1:11] <- c("Sample_ID", "Sample_ID_time", "Sample_organism", "Patient_ID", "Original_diagnosis", "disease_status", 
#                         "biopsy_location", "time", "age", "sex", "ethnicity")
# metadata_curated <- metadata[, 1:11]
# write.table(metadata_curated, '/project/gca/yuzhao1/work/final_RC2rna/bulk/metadata/metadatacurated.csv', sep = ',', row.names = F)

################ read metadata ##############
metadata <- read.table('/project/gca/yuzhao1/work/final_RC2rna/bulk/metadata/metadatacurated.csv', sep = ',', header = T)
metadata$Patient_ID <- gsub('patient: ', '', metadata$Patient_ID)
metadata$Original_diagnosis <- gsub('diagnosis: ', '', metadata$Original_diagnosis)
metadata$disease_status <- gsub('prognosis: ', '', metadata$disease_status)
metadata$biopsy_location <- gsub('tissue: ', '', metadata$biopsy_location)
metadata$time <- gsub('biopsytime: ', '', metadata$time)
metadata$age <- gsub('age: ', '', metadata$age)
metadata$sex <- gsub('Sex: ', '', metadata$sex)
metadata$ethnicity <- gsub('ethnicity: ', '', metadata$ethnicity)
rownames(metadata) <- metadata$Sample_ID


































