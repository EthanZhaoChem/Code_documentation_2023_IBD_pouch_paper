dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(plyr)
library(dplyr)
library(stringr)
source('/project/gca/yuzhao1/work/manu/rc2/scripts/rna_annotation_markers.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')

# 1. read data
out.dir <- '~/yuzhao1/work/final_RC2rna/0revision/reviewer3/4TF/'
proj <- loadArchRProject(path = "~/yuzhao1/work/final_RC2atac/dataset/gut_2kmin_6TSS_DoubletRatio2_epithelial_filtered2_EC/")
motifPositions <- getPositions(proj)

# 2, select motifs
motifs <- c('CDX1', 'FOXP1', 'KLF5', 'EHF',  'NFIA', 'GATA6', 'BACH1',
            'ESRRG', 'MAF', 'TBX3', 'HNF4G',  'NR1H4', 'PPARA')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% c("MAFA_146", "MAFF_147", "MAFG_148", "MAFK_149", "MAFB_150")]
markerMotifs

# calculate footprints
proj <- addGroupCoverages(ArchRProj = proj, 
                          groupBy = "anno1_loc",
                          useLabels = T,
                          minReplicates = 12,
                          maxReplicates = 24,
                          minCells = 40,
                          maxCells = 10000,
                          force = T)
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "anno1_loc"
)

# plot foot prints
custom_colors_EC_location <- c(
  'EC-AC'='#35978f', 'EC-POU2'='#fdb462', 'EC-TI'='#4292c6', 'EC-PP'='#9ecae1', 'EC-POU1'='#b3de69'
)
plotFootprints(
  seFoot = seFoot,
  pal = custom_colors_EC_location,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  pal = custom_colors_EC_location,
  ArchRProj = proj, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plot.name <- 'Footprints-Subtract-Bias.pdf'
file.copy(from = paste0(proj@projectMetadata$outputDirectory, '/Plots/', plot.name),
          to = paste0(out.dir, plot.name), overwrite = T)
plot.name <- 'Footprints-Divide-Bias.pdf'
file.copy(from = paste0(proj@projectMetadata$outputDirectory, '/Plots/', plot.name),
          to = paste0(out.dir, plot.name), overwrite = T)







