dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(purrr)
library(stringr)
library(harmony)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')

# path
folder.name.out <- '~/yuzhao1/work/final_RC2rna/preprocessing/4QC_plots/'



# 1. read raw file
gut_raw <- readRDS('~/gca/yuzhao1/work/final_RC2rna/preprocessing/RC2rna_all24samples_removedAmbientRNA_calculatedDoubletScores_seurat.rds')

# 2. check quality
# calculate MT and RP percentage
gut_raw[["percent.mt"]] <- PercentageFeatureSet(gut_raw, pattern = "^MT-")
gut_raw[["percent.rb"]] <- PercentageFeatureSet(gut_raw, pattern = "^RP[SL]")

# log10GenesPerUMI won't be used because plasma cell express lots of IG genes (usually require higher than 0.8)
gut_raw$log10GenesPerUMI <- log10(gut_raw$nFeature_RNA) / log10(gut_raw$nCount_RNA) 

# check the distribution stats
gut_raw@meta.data[, c("nCount_RNA", "nFeature_RNA")] %>% map(log10) %>% 
  map(~c(10^(mean(.x) + 3*sd(.x)), 10^(mean(.x) - 3*sd(.x))))%>%
  print()

# plot sample-wise QC
sel <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "log10GenesPerUMI", "Doublet_score")

for (feature in sel){
  png(paste0(folder.name.out , 'raw_', feature, '.png'),res = 300, height = 1200, width = 1800)
  df <- gut_raw@meta.data
  # if the feature is a gene, should retrive the data and add to df
  if(!(feature %in% colnames(gut_raw@meta.data))){
    addinfo <- FetchData(gut_raw, feature)
    df <- bind_cols(gut_raw@meta.data, addinfo)
  }
  
  pp <- ggviolin(df, x = 'Sample_ID', y = feature,   fill = 'biopsy_location', palette = manual_colors_rc2_location,
                 scale = 'width', width=0.8, 
                 trim = T, position = position_dodge(0.7), draw_quantiles = 0.5)+
    labs( x = NULL) +
    theme_pubr()+
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10, angle=60, hjust = 1, vjust = 1),
          axis.title = element_text(size=10),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.position = "bottom",
          plot.title = element_text(size=15, hjust=0.5, face = 'bold'))+
    labs(title = feature)
  print(pp)
  dev.off()
}

# 3. set filtering parameters for cell
nFeature_RNA_min <- 200
nFeature_RNA_max <- 6000
nCount_RNA_min <- 500
nCount_RNA_max <- 20000
percent.mt.max <- 50
percent.rb.max <- 40
log10GenesPerUMI.min <- 0.7
Doublet_score.max <- 0.2

sum(gut_raw$nCount_RNA < 20000)/ncol(gut_raw)

gut_raw_filtered<-
  subset(
    gut_raw,
    subset = 
      nFeature_RNA > nFeature_RNA_min 
    & nFeature_RNA < nFeature_RNA_max
    & nCount_RNA > nCount_RNA_min
    & nCount_RNA < nCount_RNA_max
    & percent.mt < percent.mt.max
    & percent.rb < percent.rb.max
    & log10GenesPerUMI > log10GenesPerUMI.min
    & Doublet_score < Doublet_score.max
  )

# plot sample-wise QC
sel <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "log10GenesPerUMI", "Doublet_score")
df <- gut_raw_filtered@meta.data

for (feature in sel){
  png(paste0(folder.name.out , 'filtered_', feature, '.png'),res = 300, height = 1200, width = 1800)
  # if the feature is a gene, should retrive the data and add to df
  if(!(feature %in% colnames(gut_raw_filtered@meta.data))){
    addinfo <- FetchData(gut_raw_filtered, feature)
    df <- bind_cols(gut_raw_filtered@meta.data, addinfo)
  }
  
  pp <- ggviolin(df, x = 'Sample_ID', y = feature,   fill = 'biopsy_location', palette = manual_colors_rc2_location,
                 scale = 'width', width=0.8, 
                 trim = T, position = position_dodge(0.7), draw_quantiles = 0.5)+
    labs( x = NULL) +
    theme_pubr()+
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10, angle=60, hjust = 1, vjust = 1),
          axis.title = element_text(size=10),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.position = "bottom",
          plot.title = element_text(size=15, hjust=0.5, face = 'bold'))+
    labs(title = feature)
  print(pp)
  dev.off()
}


# 4 filter genes

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = gut_raw_filtered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 120 TRUE values per gene (average 5/sample)
keep_genes <- Matrix::rowSums(nonzero) >= 120

sum(keep_genes)
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
final_filtered <- CreateSeuratObject(filtered_counts, meta.data = gut_raw_filtered@meta.data)

# 5. save filtered file
saveRDS(final_filtered, '~/yuzhao1/work/final_RC2rna/preprocessing/RC2rna_all24samples_removedAmbientRNA_calculatedDoubletScores_seurat_filtered.rds')








