dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/helper_archr.R')
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/rna_deg_markers.R')
source('~/yuzhao1/work/manu/rc2/scripts/tfs.R')

###############################################################################
# 1. preparation
addArchRThreads(4)
tf_names <- readRDS('~/yuzhao1/work/manu/rc2/plots/6tf_logo/tf_names.rds')
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2/")
union <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
seurat <- subset(union, cells = Cells(union)[union$anno1 %in% c('EC1-1', 'EC1-2', "EC2-1", "EC2-2")])
df_customized_enrichment <- readRDS("/project/gca/yuzhao1/work/final_RC2atac/peaks/customized_enrichment/EC_allPossible_enrichments.rds")

library(openxlsx)  
openxlsx::write.xlsx(df_customized_enrichment, file = paste0('~/yuzhao1/work/manu/rc2/tables/TF_EC_allPossible_enrichments.xlsx'), overwrite = T)

out.dir <- '/project/gca/yuzhao1/work/final_RC2atac/peaks/7motif_finalization/'


seurat.list <- list()
for (temp.anno in unique(seurat$anno1_loc)){
  seurat.list[[temp.anno]] <- subset(seurat, anno1_loc == temp.anno)
}

genes.broadlyExpressed.pool <- list()
for (temp.anno in unique(seurat$anno1_loc)){
  seurat.temp <- seurat.list[[temp.anno]]
  ncellsPerGene <- rowSums2(seurat.temp@assays$RNA@counts > 0.5) # cutoff is 0.5 gene for a cell (consider >0.5 as 1)
  genes.broadlyExpressed.idx <- order(ncellsPerGene/ncol(seurat.temp), decreasing = T)[1:3000] # top5000
  temp <- rownames(seurat.temp)[genes.broadlyExpressed.idx]
  genes.broadlyExpressed.pool[[temp.anno]] <- temp
}
genes.broadlyExpressed.pool_top3k <- unique(unlist(genes.broadlyExpressed.pool))


###############################################################################
# 1. selected top TFs (mlog10adjP>3, top 30)
selected_contrasts <- c("pouch2vsAC_excludingTIvsAC",
                        "pouch2_only",  # POU2vsTI minus ACvsTI
                        "EC_PPvsPOU2",
                        "EC_POU2vsPP",
                        "EC_ACvsTI",
                        "EC_TIvsAC",
                        'colon_core',
                        'ileum_core',
                        'EC_POU2vsAC',
                        'EC_ACvsPOU2',
                        'EC_POU2vsPOU1',
                        'EC_POU1vsPOU2'
)

selected_tfs_list <- list()
for (single.contrast in selected_contrasts){
  df <- df_customized_enrichment[[single.contrast]]
  temp1 <-   df[df$mlog10Padj > 3, ]
  temp1 <- temp1[head(order(temp1$mlog10Padj, decreasing = T), 1000), ]
  selected_tfs_list[[single.contrast]] <- temp1$feature
}

selected_tfs <- unique(unlist(selected_tfs_list))
selected_tfs <- selected_tfs %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,1)


# 2. filter for broadly expressed genes
tfs <- intersect(selected_tfs, genes.broadlyExpressed.pool_top3k)
cat(sort(tfs))

# 3. filter for positively regulating TFs (>0.1), this gives a broad list of genes
positive_regulators_EC_corr0.1 <- readRDS('~/yuzhao1/work/manu/rc2/plots/7atac_positiveTFregulator/positive_regulators_EC_corr0.1.rds')
tfs_positive <- intersect(tfs, positive_regulators_EC_corr0.1)

negative_regulators_EC_corr0.5 <- readRDS('~/yuzhao1/work/manu/rc2/plots/7atac_negativeTFregulator/negative_regulators_EC_corr0.5.rds')
tfs_negative <- intersect(tfs, negative_regulators_EC_corr0.5)

sort(tfs_positive)
tfs_final <- tfs_positive




# 4. use the finalized TFs to build the df
df <- tf_names
for(single.contrast in selected_contrasts){
  df_stat <- df_customized_enrichment[[single.contrast]]
  df[[single.contrast]] <- df_stat[df$tf_full, 'mlog10Padj']
}

tfs_final <- c(tf_genes_negative, tf_genes_refined)
# set the non-significant values to zero (induced by TFs in other celltypes)
df_filtered <- df[df$tf%in%tfs_final,]
for (i in 1:nrow(df_filtered)) {
  for(j in 3:ncol(df_filtered)){
    df_filtered[i,j] <- ifelse(df_filtered[i,j]<3, 0, df_filtered[i,j])
  }
}
View(df_filtered)

saveRDS(df_filtered, paste0(out.dir, 'df_finalized_enriched_motifs.rds'))
saveRDS(df_filtered, paste0(out.dir, 'df_finalized_enriched_motifs_negative.rds')) # saved the negatively correlated TFs

# then in manu folder, plot heatmap of gene expression, chromVAR
# plot vln plot of DEGs (select by observation)
# plot chromVAR score for heterogeneity
# pseudo time, logo, in silico test for selected motifs
# negative/positive regulator characterization
# use DORC to build a network (see which are DORCs among these selected TFs), can i use cell oracle to build network?






