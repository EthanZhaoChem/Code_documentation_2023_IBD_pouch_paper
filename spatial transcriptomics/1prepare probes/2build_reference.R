dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(harmony)
library(Seurat)
library(biomaRt)
library(org.Hs.eg.db)
source('~/yuzhao1/scripts/plot.R')
out.dir <- '~/yuzhao1/work/final_RC2rna/0revision/spatial/ref'

##################  ################ ##################  ################
# 1. read data
gut_all <- readRDS('~/yuzhao1/work/final_RC2rna/annotation/rds/gut_filtered.rds')
cells_pouch <- Cells(gut_all)[which(gut_all$biopsy_location =='POU')]
seurat_obj <- subset(gut_all, cells = cells_pouch)

# 2, round to integers
seurat_obj@assays$RNA@counts <- round(seurat_obj@assays$RNA@counts)
all.equal(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x, as.integer(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x))
head(GetAssayData(object=seurat_obj, assay="RNA", slot="counts")@x)

# 3, subsample if necessary (max data size should be 500MB)
subsample_rate <- 1
subset <- 1:ncol(seurat_obj)
if (subsample_rate < 1) {
  subset <- sample(subset, subsample_rate * length(subset))
}

# 4, add Gene ID
gene_symbols <- seurat_obj@assays$RNA@meta.features %>% rownames()
gdf <- biomaRt::select(org.Hs.eg.db, gene_symbols, "ENSEMBL", "SYMBOL")
gdf <- gdf[!base::duplicated(gdf$SYMBOL),]
seurat_obj@assays$RNA@meta.features$Gene_symbol <- gene_symbols
seurat_obj@assays$RNA@meta.features$Gene_ensembl <- mapvalues(gene_symbols, gdf$SYMBOL, gdf$ENSEMBL, warn_missing = F)
seurat_obj@assays$RNA@meta.features$feature_type <- 'Gene Expression'

# 5, write data to file
writeCounts <- function(out_dir, data, barcodes = colnames(data), gene.id = rownames(data), gene.symbol = rownames(data), feature.type = "Gene Expression", subset = 1:length(barcodes)) {
  require("R.utils")
  require("Matrix")
  
  if (file.exists(out_dir) || (dir.exists(out_dir) && length(list.files(out_dir)) > 0)) {
    stop("The specified output directory already exists! Not overwriting")
  }
  dir.create(out_dir, recursive = TRUE)
  
  if (require("data.table")) {
    data.table::fwrite(
      data.table::data.table(barcode = barcodes[subset]),
      file.path(out_dir, "barcodes.tsv.gz"),
      col.names = FALSE
    )
    
    data.table::fwrite(
      data.table::data.table(
        gene_id = gene.id,
        gene_symbol = gene.symbol,
        feature_type = feature.type
      ),
      file.path(out_dir, "features.tsv.gz"),
      col.names = FALSE
    )
  } else {
    write.table(
      data.frame(barcode = barcodes[subset]),
      gzfile(file.path(out_dir, "barcodes.tsv.gz")),
      sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE
    )
    
    write.table(
      data.frame(
        gene_id = gene.id,
        gene_symbol = gene.symbol,
        feature_type = feature.type
      ),
      file.path(out_dir, "features.tsv.gz"),
      sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE
    )
  }
  
  Matrix::writeMM(data[, subset], file.path(out_dir, "matrix.mtx"))
  R.utils::gzip(file.path(out_dir, "matrix.mtx"), remove = TRUE)
}

# Run function
writeCounts(
  out_dir = out.dir,
  GetAssayData(seurat_obj, assay="RNA", slot="counts"),
  gene.id = seurat_obj@assays$RNA@meta.features$Gene_ensembl,
  gene.symbol = seurat_obj@assays$RNA@meta.features$Gene_symbol,
  feature.type = GetAssay(seurat_obj)@meta.features[["feature_type"]],
  barcodes = colnames(seurat_obj), 
  subset = subset
  )
list.files("~/yuzhao1/work/final_RC2rna/0revision/spatial/ref")

# 6. save annotation
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
  if (require("data.table", quietly = TRUE)) {
    data.table::fwrite(
      data.table::data.table(
        barcode = barcodes,
        annotation = unlist(seurat_obj[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv")
    )
  } else {
    write.table(
      data.frame(
        barcode = barcodes,
        annotation = unlist(seurat_obj[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv"),
      sep = ",", row.names = FALSE
    )
  }
  
  bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))
  
  utils::zip(
    bundle,
    list.files(out_dir, full.names = TRUE),
    zip = "zip"
  )
  
  if (file.info(bundle)$size / 1e6 > 500) {
    warning("The output file is more than 500 MB and will need to be subset further.")
  }
}

# Run function
bundleOutputs(out_dir = out.dir, data = seurat_obj, subset = subset, cell_type = 'anno2')
list.files("~/yuzhao1/work/final_RC2rna/0revision/spatial/ref")






