library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('/project/gca/yuzhao1/work/manu/rc2/scripts/rna_annotation_markers.R')
source('/project/gca/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')

# 1. old version
tfs <- read.table('~/yuzhao1/resource/scenic/utoronto_human_tfs_v_1.01.txt', header = F)[[1]]

velogenes.all <- read.table('~/yuzhao1/work/final_RC2rna/velocity/epithelial/velocity_genes.csv', header = F)[[1]]
velogenes.ec_stem <- read.table('~/yuzhao1/work/final_RC2rna/velocity/epithelial/ec_stem/velocity_genes.csv', header = F)[[1]]
velogenes.pouch_ec_stem <- read.table('~/yuzhao1/work/final_RC2rna/velocity/epithelial/pouch_ec_stem/velocity_genes.csv', header = F)[[1]]

intersect(tfs, velogenes.all) %>% sort()
intersect(tfs, velogenes.ec_stem) %>% sort()
intersect(tfs, velogenes.pouch_ec_stem) %>% sort()

# run EC markers in ~/yuzhao1/work/pouch_atac2/R/Figures/1Celltype.R
c(velogenes.all, velogenes.ec_stem, velogenes.pouch_ec_stem) %>% unique() %>% intersect(., EC_markers)
# get ("MEIS1"   "TFCP2L1" "SATB2"   "PPARG"   "ZBTB10"  "PAX6"    "CREB3L1" "SMAD9"   "NR1H4" )

# 2. new version
velogenes.anno1 <- read.table('~/yuzhao1/work/final_RC2rna/velocity/epithelial/velocity_genes_differential_by_anno1.csv', header = T, sep = ',')
c(velogenes.anno1$EC1.1, velogenes.anno1$EC2.1) %>% intersect(., rna_deg_markers_ec_tfs)


pouch_velogenes.anno1 <- read.table('~/yuzhao1/work/final_RC2rna/velocity/epithelial/pouch_velocity_genes_differential_by_anno1.csv', header = T, sep = ',')
c(pouch_velogenes.anno1$EC1.1, pouch_velogenes.anno1$EC2.1) %>% intersect(., rna_deg_markers_ec_tfs)
