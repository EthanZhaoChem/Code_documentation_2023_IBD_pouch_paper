dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(Seurat)
library(Matrix)
library(magrittr)
SampleIDs <- read.csv("~/yuzhao1/work/final_RC2rna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID

mtx_outdir <- '~/yuzhao1/work/final_RC2rna/preprocessing/matrix_removedAmbientRNA/'
gene_outdir <- '~/yuzhao1/work/final_RC2rna/preprocessing/gene_removedAmbientRNA/'
barcodes_outdir <- '~/yuzhao1/work/final_RC2rna/preprocessing/barcodes_outdir/'
doublet_outdir <- '~/yuzhao1/work/final_RC2rna/preprocessing/doublet_scores/'

mtx_outdir.gca <- '~/yuzhao1/work/final_GCArna/preprocessing/matrix_removedAmbientRNA/'
gene_outdir.gca <- '~/yuzhao1/work/final_GCArna/preprocessing/gene_removedAmbientRNA/'
barcodes_outdir.gca <- '~/yuzhao1/work/final_GCArna/preprocessing/barcodes_outdir/'
doublet_outdir.gca <- '~/yuzhao1/work/final_GCArna/preprocessing/doublet_scores/'
out.srt.obj.path <- '~/yuzhao1/work/final_RC2rna/preprocessing/RC2rna_all24samples_removedAmbientRNA_calculatedDoubletScores_seurat.rds'

# read each sample matrix (24 samples for RC2 analysis) and combine as a seurat object
srat = list()
doubletscore.all <- c()
sample.col <- c()
for (i in 1:24) {
  SampleID <- SampleIDs[[i]]
  
  # read samples in RC2
  if(grepl('PP', SampleID)|grepl('POU', SampleID)){
    mtx <- readMM(paste0(mtx_outdir, SampleID,'.mtx'))
    genes <- read.table(file = paste0(gene_outdir, SampleID,'.tsv'),
                        sep = '\t', header = F) %>% unlist()
    cells <- read.table(file = paste0(barcodes_outdir, SampleID,'.tsv'),
                        sep = '\t', header = F) %>% unlist()
    doublet.scores <- read.table(file = paste0(doublet_outdir, SampleID,'.csv'),
                                 sep = ',', header = F) %>% unlist()
    names(doublet.scores) <- cells
    colnames(mtx) <- cells
    rownames(mtx) <- genes
    
    doubletscore.all <- c(doubletscore.all, doublet.scores)
    sample.col <- c(sample.col, rep(SampleID, length(cells)))
    srat[[SampleID]] <- mtx
  }
  
  # read samples in GCA
  if(grepl('TI', SampleID)|grepl('AC', SampleID)){
    mtx <- readMM(paste0(mtx_outdir.gca, SampleID,'.mtx'))
    genes <- read.table(file = paste0(gene_outdir.gca, SampleID,'.tsv'),
                        sep = '\t', header = F) %>% unlist()
    cells <- read.table(file = paste0(barcodes_outdir.gca, SampleID,'.tsv'),
                        sep = '\t', header = F) %>% unlist()
    doublet.scores <- read.table(file = paste0(doublet_outdir.gca, SampleID,'.csv'),
                                 sep = ',', header = F) %>% unlist()
    names(doublet.scores) <- cells
    colnames(mtx) <- cells
    rownames(mtx) <- genes
    
    doubletscore.all <- c(doubletscore.all, doublet.scores)
    sample.col <- c(sample.col, rep(SampleID, length(cells)))
    srat[[SampleID ]] <- mtx
  }

}


# Combine all count matricies into one matrix
srat.mtx = do.call(cbind, srat)
srat.obj = CreateSeuratObject(srat.mtx)



# add metadata from processing to seurat object
srat.obj$Sample_ID <- sample.col
srat.obj$Doublet_score <- doubletscore.all
srat.obj$Patient_ID <- 
  srat.obj$Sample_ID %>%
  gsub('-', '', .) %>%
  gsub('AC', '', .) %>%
  gsub('TI', '', .) %>%
  gsub('PP', '', .) %>%
  gsub('POU', '', .) 

srat.obj$biopsy_location <- 'na'

srat.obj$biopsy_location[grep('AC', srat.obj$Sample_ID)] <- 'AC'
srat.obj$biopsy_location[grep('TI', srat.obj$Sample_ID)] <- 'TI'
srat.obj$biopsy_location[grep('PP', srat.obj$Sample_ID)] <- 'PP'
srat.obj$biopsy_location[grep('POU', srat.obj$Sample_ID)] <- 'POU'

saveRDS(srat.obj, out.srt.obj.path)

















