dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(MuSiC)
library(Biobase)
library(dplyr)
library(plyr)
library(reshape2)
library(cowplot)
library(SingleCellExperiment)
library(Seurat)
setwd('/project/gca/yuzhao1/work/final_RC2rna/bulk/music/dataset')

# in this script, single end libs are counted by reads
# paired end libs are counted by fragements



################ read reference data ##############
ref <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
cells_pppou <- Cells(ref)[which(ref$biopsy_location %in% c('PP', 'POU'))]
ref_pppou <- subset(ref, cells = cells_pppou)

# cells_pou <- Cells(ref)[which(ref$biopsy_location %in% c('POU'))]
# ref_pou <- subset(ref, cells = cells_pou)

################ reference build up ##############
ref_music<- ref_pppou
# xx <- unique(ref_music$anno1)
# for (i in 1:length(xx)) {
#   paste0("'", xx[[i]], "'", " = '', \n") %>% cat()
# }

labels_music <- list(
    'Stem' = 'Progenitor', 
    'TA' = 'Progenitor', 
    'M-like' = 'Rare cells', 
    'EC1-1' = 'EC1', 
    'EC1-2' = 'EC1', 
    'EC2-1' = 'EC2', 
    'EC2-2' = 'EC2', 
    'Goblet1' = 'Goblet', 
    'Goblet2' = 'Goblet', 
    'EEC' = 'Rare cells', 
    'Paneth' = 'Rare cells', 
    'Tuft' = 'Rare cells', 
    'BEST4' = 'Rare cells', 
    'CD4 Tcm' = 'Lymphocyte', 
    'CD103- CD4 Trm' = 'Lymphocyte', 
    'CD103+ CD4 Trm' = 'Lymphocyte', 
    'Treg' = 'Lymphocyte', 
    'CD103+ CD8 Trm' = 'Lymphocyte', 
    'KLRG1+ CD8 Trm' = 'Lymphocyte', 
    'MAIT' = 'Lymphocyte', 
    'NK' = 'Lymphocyte', 
    'NK T' = 'Lymphocyte', 
    'ILCs' = 'Lymphocyte', 
    'gdT' = 'Lymphocyte', 
    'IgA plasma' = 'Plasma', 
    'IgG plasma' = 'Plasma', 
    'Naive B' = 'Lymphocyte', 
    'Memory B' = 'Lymphocyte', 
    'GC B' = 'Lymphocyte', 
    'Monocyte' = 'Myeloid', 
    'Macrophage' = 'Myeloid', 
    'Neutrophil' = 'Myeloid', 
    'Mast' = 'Myeloid', 
    'cDC2' = 'Myeloid', 
    'Lymphoid DC' = 'Myeloid', 
    'cDC1' = 'Myeloid', 
    'Venous' = 'Stromal', 
    'Stromal-1' = 'Stromal', 
    'Arterial' = 'Stromal', 
    'Stromal-2' = 'Stromal', 
    'Lymphatic endothelium' = 'Stromal', 
    'Pericyte' = 'Stromal', 
    'Myofibroblast' = 'Stromal', 
    'Glial' = 'Stromal', 
    'Stromal-3' = 'Stromal', 
    'Contractile pericyte' = 'Stromal', 
    'Smooth muscle' = 'Stromal')
ref_music$anno_music <- mapvalues(ref_music$anno1, names(labels_music), labels_music) %>% unlist()

################ read metadata ##############
metadata <- read.table('/project/gca/yuzhao1/work/final_RC2rna/bulk/metadata/metadatacurated.csv', sep = ',', header = T)
metadata$SRRID <- 3493774:3493850
metadata$Seq_Method <- 'Paired End'
metadata[which(metadata$SRRID < 3493790), 'Seq_Method'] <- 'Single End'

metadata$Patient_ID <- gsub('patient: ', '', metadata$Patient_ID)
metadata$Original_diagnosis <- gsub('diagnosis: ', '', metadata$Original_diagnosis)
metadata$disease_status <- gsub('prognosis: ', '', metadata$disease_status)
metadata$biopsy_location <- gsub('tissue: ', '', metadata$biopsy_location)
metadata$time <- gsub('biopsytime: ', '', metadata$time)
metadata$age <- gsub('age: ', '', metadata$age)
metadata$sex <- gsub('Sex: ', '', metadata$sex)
metadata$ethnicity <- gsub('ethnicity: ', '', metadata$ethnicity)
rownames(metadata) <- metadata$Sample_ID


################ prepare bulk counts ##############
counts <- readRDS('/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/5Final_counts_matrix.rds')
colnames(counts) <- metadata$Sample_ID # raw metadata's sample order is same with the SRR order

################ workflow ##############
# 1. bulk data
metadata_description <- data.frame(labelDescription= c("Sample_ID", "Sample_ID_time", "Sample_organism", "Patient_ID", 
                                                       "Original_diagnosis", "disease_status", "biopsy_location", "time", 
                                                       "age", "sex", "ethnicity"), 
                                   row.names=c("Sample_ID", "Sample_ID_time", "Sample_organism", "Patient_ID", 
                                               "Original_diagnosis", "disease_status", "biopsy_location", "time", 
                                               "age", "sex", "ethnicity"))
SC.eset = ExpressionSet(assayData = data.matrix(counts))
bulk.mtx = exprs(SC.eset)

# 2. single cell data
gene_annotations <- data.frame(GeneID = rownames(ref_music),
                               GeneName = rownames(ref_music))
rownames(gene_annotations) <- rownames(ref_music)

matrix.sce <- as.matrix(ref_music@assays$RNA@counts) 
rownames(matrix.sce) <- NULL
colnames(matrix.sce) <- NULL

ref.sce <- SingleCellExperiment(assays= list(counts = matrix.sce),
                                colData = ref_music@meta.data,
                                rowData = gene_annotations)

# 3. Estimate cell type proportions
Est.prop = music_prop(
  bulk.mtx = bulk.mtx,
  sc.sce = ref.sce,
  clusters = 'anno_music',
  samples = 'Sample_ID',
  select.ct = unique(ref_music$anno_music),
  cell_size = NULL, # estimate cell size from data
  centered = F,
  normalize = F,
  verbose = T
)
names(Est.prop)

saveRDS(Est.prop, '/project/gca/yuzhao1/work/final_RC2rna/bulk/music/Est.prop_pppou.rds')


Est.prop <- readRDS('/project/gca/yuzhao1/work/final_RC2rna/bulk/music/Est.prop_pppou.rds')

results <- merge(Est.prop$Est.prop.weighted, metadata, by = "row.names")
colnames(results)[[1]] <- 'Sample_ID'
rownames(results) <- results$Sample_ID
write.csv(results, '/project/gca/yuzhao1/work/final_RC2rna/bulk/music/bulk_deconvolution_withPPandPOU_asRef.csv')
saveRDS(results, '/project/gca/yuzhao1/work/final_RC2rna/bulk/music/Final_results_bulk_deconvolution_withPPandPOU_asRef.rds')

# 
# a1 <- xx[, 'EC2'] %>% sort(decreasing = T) %>% head(20) %>% names()
# a2 <- metadata$Sample_ID[which(metadata$biopsy_location == 'Ileal pouch')]
# intersect(a1,a2) %>% length()









































