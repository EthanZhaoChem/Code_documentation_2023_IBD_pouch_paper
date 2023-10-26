# 12 quantitative colors
library(ggsci)
library(ggpubr)
library(colorspace)
library(paletteer) 
library(wesanderson)
library(ggrepel)
library(ggrastr)
library(ggforce)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(scales)



############ 1, rna umap ##############
custom_colors_rna_lineage.levels <- c("epithelial", "tcell", "bcell", "plasma cell", "myeloid", "others")
custom_colors_rna_epithelial.levels <- c("Stem", "TA", "EC1-1", "EC1-2", "EC2-1", "EC2-2", "Goblet1",
                                         "Goblet2", "M-like", "BEST4", "Paneth", "EEC", "Tuft")
custom_colors_rna_tcell.levels <- c("CD4 Tcm" , "Treg", "CD103- CD4 Trm", "CD103+ CD4 Trm", 
                                     "CD103+ CD8 Trm", "KLRG1+ CD8 Trm", "gdT", 
                                     "MAIT", "NK T", "NK", "ILCs" )
custom_colors_rna_bcell.levels <- c("GC B", "Naive B", "Memory B", "IgA plasma", "IgG plasma")
custom_colors_rna_myeloid.levels <- c("Monocyte", "Macrophage", "cDC1", "cDC2", "Lymphoid DC", "Mast", "Neutrophil")
custom_colors_rna_others.levels <- c("Stromal-1", "Stromal-2", "Stromal-3", "Myofibroblast", "Arterial", "Venous", "Pericyte", 
                                     "Contractile pericyte" , "Smooth muscle", "Lymphatic endothelium" , "Glial" )


custom_colors_rna_lineage <- c('epithelial' = '#72b9bf',
                               'tcell' = '#82c07f',
                               'bcell' = '#e89f9e',
                               'myeloid' = '#cbc8e1',
                               'others' = '#efb788')
custom_colors_rna_epithelial <- c("Stem" = '#fccde5', #purple
                                  "TA" = '#82c07f', # green
                                  "EC1-1" = '#9ecae1',
                                  "EC1-2" = '#4292c6',
                                  "EC2-1" = '#80cdc1',
                                  "EC2-2" = '#35978f',
                                  "Goblet1" = '#fbb4ae',
                                  "Goblet2" = '#fed9a6',
                                  "M-like" = '#8073ac',
                                  "BEST4" = '#7998bb',
                                  "Paneth" = '#cd62af',
                                  "EEC" = '#fb8072',
                                  "Tuft" = '#e6f598')

custom_colors_rna_epithelial_anno2 <- c("Stem1" = '#fccde5', 
                                        "Stem2" = '#e7298a',
                                        "TA" = '#82c07f', # green
                                        "EC1-1" = '#9ecae1',
                                        "EC1-2" = '#4292c6',
                                        "EC2-1" = '#80cdc1',
                                        "EC2-2" = '#35978f',
                                        "Goblet1" = '#fbb4ae',
                                        "Goblet2" = '#fed9a6',
                                        "M-like" = '#8073ac',
                                        "BEST4" = '#7998bb',
                                        "Paneth" = '#cd62af',
                                        "EEC" = '#fb8072',
                                        "Tuft" = '#e6f598')

custom_colors_rna_tcell <- c("CD4 Tcm" = '#ffed6f', # yellow
                             "Treg" = '#bebada', # light purple
                             "CD103- CD4 Trm" = '#b2df8a', # green
                             "CD103+ CD4 Trm" = '#33a02c', # green
                             "CD103+ CD8 Trm" = '#fc8d59', # orange
                             "KLRG1+ CD8 Trm" = '#fc4e2a', # orange
                             "gdT" = '#bc80bd', # purple
                             "MAIT" = '#4393c3', # blue
                             "NK T" = '#c51b7d', # red purple
                             "NK" = '#de77ae', # red purple
                             "ILCs" = '#8c510a') # brown


custom_colors_rna_bcell <- c("GC B" = '#ccebc5', # blue green
                             "Naive B" = '#b3de69', # green
                             "Memory B" = '#80b1d3', # blue
                             "IgA plasma" = '#8dd3c7', # blue green
                             "IgG plasma" = '#fb8072') # # cute red


custom_colors_rna_myeloid <- c("Monocyte" = '#6fbe6d', # green
                             "Macrophage" = '#4393c3', # blue
                             "cDC1" = '#542788', # purple
                             "cDC2" = '#d3539d', # purple
                             "Lymphoid DC" = '#01665e', # # purple
                             "Mast" = '#cfa061', # brown
                             "Neutrophil" = '#db433a') # cute red


custom_colors_rna_others <- c("Stromal-1" = '#8dd3c7',
                              "Stromal-2" = '#fb8072',
                              "Stromal-3" = '#bebada',
                              "Myofibroblast" = '#80b1d3',
                              "Arterial" = '#fdb462',
                              "Venous" = '#b3de69',
                              "Pericyte" = '#fccde5',
                              "Contractile pericyte" = '#bc80bd',
                              "Smooth muscle" = '#ccebc5',
                              "Lymphatic endothelium" = '#ffed6f',
                              "Glial" = '#ec5f60')


custom_colors_rna_immune <- c(custom_colors_rna_tcell, custom_colors_rna_bcell, custom_colors_rna_myeloid)

custom_colors_rna_all <- c(custom_colors_rna_tcell, custom_colors_rna_bcell, custom_colors_rna_myeloid,
                           custom_colors_rna_epithelial, custom_colors_rna_others)


############ 2, atac umap ##############
custom_colors_atac_lineage.levels <- c("epithelial", "tcell", "bcell", "plasma cell", "myeloid", "others")
custom_colors_atac_epithelial.levels <- c("Stem", "TA", "EC1-1", "EC1-2", "EC2-1", "EC2-2", "Goblet1",
                                         "Goblet2", "M-like", "BEST4", "Paneth", "EEC", "Tuft")
custom_colors_atac_tcell.levels <- c("CD4 Tcm" , "Treg", "CD103- CD4 Trm", "CD103+ CD4 Trm", 
                                    "CD103+ CD8 Trm", "KLRG1+ CD8 Trm", "gdT", 
                                    "MAIT", "NK T", "NK", "ILCs" )
custom_colors_atac_bcell.levels <- c("GC B", "Naive B", "Memory B", "Plasma cell")
custom_colors_atac_myeloid.levels <- c("Monocyte", "Macrophage", "DC", "Lymphoid DC", "Mast", "Neutrophil")
custom_colors_atac_others.levels <- c("Stromal-1", "Stromal-2", "Stromal-3", "Myofibroblast", "Arterial", "Venous", "Pericyte", 
                                     "Contractile pericyte" , "Smooth muscle", "Lymphatic endothelium" , "Glial" )

custom_colors_atac_lineage <- c('epithelial' = '#72b9bf',
                               'tcell' = '#82c07f',
                               'bcell' = '#e89f9e',
                               'myeloid' = '#cbc8e1',
                               'others' = '#efb788')

custom_colors_atac_epithelial <- c("Stem" = '#fccde5', #purple
                                  "EC1-1" = '#9ecae1',
                                  "EC1-2" = '#4292c6',
                                  "EC2-1" = '#80cdc1',
                                  "EC2-2" = '#35978f',
                                  "Goblet1" = '#fbb4ae',
                                  "Goblet2" = '#fed9a6',
                                  "BEST4" = '#7998bb',
                                  "Paneth" = '#cd62af',
                                  "EEC" = '#fb8072',
                                  "Tuft" = '#e6f598')

custom_colors_atac_epithelial_anno2 <- c("Stem1" = '#fccde5', 
                                   "Stem2" = '#e7298a',
                                   "EC1-1" = '#9ecae1',
                                   "EC1-2" = '#4292c6',
                                   "EC2-1" = '#80cdc1',
                                   "EC2-2" = '#35978f',
                                   "Goblet1" = '#fbb4ae',
                                   "Goblet2" = '#fed9a6',
                                   "BEST4" = '#7998bb',
                                   "Paneth" = '#cd62af',
                                   "EEC" = '#fb8072',
                                   "Tuft" = '#e6f598')



custom_colors_atac_tcell <- c("CD4 Tcm" = '#ffed6f', # yellow
                             "Treg" = '#bebada', # light purple
                             "CD103- CD4 Trm" = '#b2df8a', # green
                             "CD103+ CD4 Trm" = '#33a02c', # green
                             "CD103+ CD8 Trm" = '#fc8d59', # orange
                             "KLRG1+ CD8 Trm" = '#fc4e2a', # orange
                             "gdT" = '#bc80bd', # purple
                             "MAIT" = '#4393c3', # blue
                             "NK" = '#de77ae', # red purple
                             "ILCs" = '#8c510a') # brown


custom_colors_atac_bcell <- c("GC B" = '#ccebc5', # blue green
                             "Naive B" = '#b3de69', # green
                             "Memory B" = '#80b1d3', # blue
                             "Plasma cell" = '#8dd3c7') # blue green


custom_colors_atac_myeloid <- c("Monocyte" = '#6fbe6d', # green
                               "Macrophage" = '#4393c3', # blue
                               "DC" = '#542788', # purple
                               "Mast" = '#cfa061', # brown
                               "Neutrophil" = '#db433a') # cute red


custom_colors_atac_others <- c("Stromal-1" = '#8dd3c7',
                              "Stromal-2" = '#fb8072',
                              "Stromal-3" = '#bebada',
                              "Myofibroblast" = '#80b1d3',
                              "Arterial" = '#fdb462',
                              "Venous" = '#b3de69',
                              "Pericyte" = '#fccde5',
                              "Lymphatic endothelium" = '#ffed6f',
                              "Glial" = '#ec5f60')


custom_colors_atac_immune <- c(custom_colors_atac_tcell, custom_colors_atac_bcell, custom_colors_atac_myeloid)

custom_colors_atac_all <- c(custom_colors_atac_tcell, custom_colors_atac_bcell, custom_colors_atac_myeloid,
                           custom_colors_atac_epithelial, custom_colors_atac_others)

############ 3, rna heat umap ##############
rc2_rna_heatmap_colors_gradient1 <-c(
  '#08306b',
  # '#08519c',
  '#2171b5',
  # '#4292c6',
  # '#6baed6',
  '#9ecae1',
  '#c6dbef',
  # '#deebf7',
  '#f7fbff',
  '#ffeda0',
  # '#fed976',
  '#feb24c',
  '#fd8d3c',
  '#fc4e2a',
  '#e31a1c',
  '#bd0026',
  '#800026'
)
# show_col(rc2_rna_heatmap_colors_gradient1)


############ 4, location color ##############
custom_colors_location <- c(
  'AC'='#80b1d3', 'TI'='#8dd3c7', 'POU'='#fdb462', 'PP'='#b3de69'
)

########### 5, atac ################

# custom_colors_atac_umap_continuous <- paletteContinuous("solarExtra")




