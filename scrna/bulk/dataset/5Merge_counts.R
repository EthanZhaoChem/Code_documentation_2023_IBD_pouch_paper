dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')

counts_matrix_singleEnd <- read.table('/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/counts_matrix_singleEnd.txt', header = T)
counts_matrix_singleEnd <- counts_matrix_singleEnd[,-c(1:6)]
colnames(counts_matrix_singleEnd) <- colnames(counts_matrix_singleEnd) %>% 
  gsub('X.project.gca.yuzhao1.work.final_RC2rna.bulk.dataset.sorted_bam.', '', .) %>%
  gsub('_sorted.bam', '', .)

counts_matrix_pairedEnd <- read.table('/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/counts_matrix_pairedEnd.txt', header = T)
counts_matrix_pairedEnd <- counts_matrix_pairedEnd[,-c(1:6)]
colnames(counts_matrix_pairedEnd) <- colnames(counts_matrix_pairedEnd) %>% 
  gsub('X.project.gca.yuzhao1.work.final_RC2rna.bulk.dataset.sorted_bam.', '', .) %>%
  gsub('_sorted.bam', '', .)

counts_matrix <-  left_join(counts_matrix_singleEnd , counts_matrix_pairedEnd , by = "gene_name")
counts_matrix$gene_name <- make.unique(counts_matrix$gene_name, '_duplicated_')
rownames(counts_matrix) <- counts_matrix$gene_name
counts_matrix$gene_name <- NULL

saveRDS(counts_matrix, '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/5Final_counts_matrix.rds')







