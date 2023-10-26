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

rna.lineage.markers <- list(
  'epithelial' = c('EPCAM', 'PHGR1','FABP1', 'KRT8'),
  'tcell' = c('CD3D', 'CCL5', 'IL7R', 'STAT4'),
  'bcell' = c('CD79A', 'MS4A1'),
  'plasma cell' = c('XBP1', 'SLAMF7', 'IGHA1', 'IGHG1'),
  'myeloid' = c('C1QA', 'LYZ', 'CPA3', 'TYROBP'),
  'others' = c('IGFBP7', 'COL3A1', 'PECAM1', 'S100B')
)


rna.epithelial.markers <- list(
  "Stem" = c('LGR5', 'OLFM4', 'SMOC2'), 
  "TA" = c('MKI67', 'TOP2A', 'TUBA1B'), 
  "EC1-1" = c('ADH1C', 'SI', 'GSTA2'),
  "EC1-2" = c('SLC15A1', 'APOA4', 'CUBN'),
  "EC2-1" = c('CD24', 'FABP5'),
  "EC2-2" = c('CEACAM5', 'CA2', 'CEACAM7', 'AQP8'),
  "Goblet1" = c('MUC2', 'TFF3', 'CLCA1', 'SPINK4', 'FER1L6'),
  "Goblet2" = c('ITLN1',  'KLK1', 'BEST2'),
  "M-like" = c('CCL20'),
  "BEST4" = c('BEST4', 'CA7', 'NOTCH2', 'SPIB'),
  "Paneth" = c('ITLN2', 'DEFA5', 'REG3A'),
  "EEC" = c('CHGA', 'KCNB2', 'RIMBP2'),
  "Tuft" = c('POU2F3', 'FYB1')
)

rna.tcell.markers <- list(
  "CD4 Tcm" = c('CD3D', 'CD4', 'LEF1', 'TCF7', 'CCR7'), 
  "Treg" = c('CTLA4', 'FOXP3', 'IL2RA'), 
  "CD103- CD4 Trm" = c('IL7R'), 
  "CD103+ CD4 Trm" = c('ITGAE', 'IL17A', 'IL26'), 
  "CD103+ CD8 Trm" = c('CD8A','CD8B'), 
  "KLRG1+ CD8 Trm" = c('KLRG1', 'GZMK','IFNG', 'HLA-DRB1', 'HLA-DRA'), 
  "gdT" = c('TRDC', 'GNLY', 'GZMA', 'ENTPD1'), 
  "MAIT" = c('SLC4A10', 'NCR3'),
  "NK T" = c('FGFBP2', 'GZMB', 'CX3CR1'),
  "NK" = c('NCR1', 'KLRF1'), 
  "ILCs" = c('PRKG1', 'PCDH9', 'AFF3', 'AREG','IL1R1', 'IL23R', 'KIT') 
)


rna.bcell.markers <- list(
  "GC B" = c('BCL6'),
  "Naive B" = c('CD19', 'MS4A1', 'IGHD'),
  "Memory B" = c('CD27', 'CD83'),
  "IgA plasma" = c('SDC1', 'SLAMF7', 'IGHA1'),
  "IgG plasma" = c('IGHG1')
)

rna.myeloid.markers <- list("Monocyte" = c('VCAN', 'FCN1', 'EREG', 'S100A8', 'S100A4'),
                            "Macrophage" = c('CD14', 'CD163', 'MMP12', 'C1QA'),
                            "cDC1" = c('CLEC9A'),
                            "cDC2" = c('CD1C'),
                            "Lymphoid DC" = c('LAMP3'),
                            "Mast" = c('CPA3', 'KIT', 'TPSB2'),
                            "Neutrophil" = c('FCGR3B')
)

rna.others.markers <- list("Stromal-1" = c('COL1A1', 'ADAMDEC1', 'CCL11', 'CCL13'), 
                           "Stromal-2" = c('NRG1', 'NPY', 'PTGS1'),  
                           "Stromal-3" = c('SOX6', 'COL4A6'), 
                           "Myofibroblast" = c('ACTA2', 'TAGLN'), 
                           "Arterial" = c('PECAM1', 'HEY1', 'EFNB2'), 
                           "Venous" = c('ACKR1', 'VWF'), 
                           "Pericyte" = c('NOTCH3', 'MCAM', 'RGS5'), 
                           "Contractile pericyte" = c('PLN', 'RERGL', 'KCNAB1'), 
                           "Smooth muscle" = c('DES', 'CNN1'), 
                           "Lymphatic endothelium" = c('PROX1', 'LYVE1', 'CCL21'), 
                           "Glial" = c('S100B', 'NRXN1'))

rna.unionHeatmap.markers <- list("Stem" = c('EPCAM','LGR5'), 
                                 "TA" = c('MKI67'), 
                                 "EC1-1" = c('GSTA1'),
                                 "EC1-2" = c('APOA4'),
                                 "EC2-1" = c('CD24'),
                                 "EC2-2" = c('CEACAM5'),
                                 "Goblet1" = c('MUC2'),
                                 "Goblet2" = c('BEST2'),
                                 "M-like" = c('CCL20'),
                                 "BEST4" = c('BEST4'),
                                 "Paneth" = c('DEFA5'),
                                 "EEC" = c('CHGA'),
                                 "Tuft" = c('FYB1'),
                                 "CD4 Tcm" = c('CD3D', 'TCF7'), 
                                 "Treg" = c('CTLA4', 'FOXP3'), 
                                 "CD103- CD4 Trm" = c('IL7R'), 
                                 "CD103+ CD4 Trm" = c('ITGAE'), 
                                 "CD103+ CD8 Trm" = c('CD8A'), 
                                 "KLRG1+ CD8 Trm" = c('KLRG1'), 
                                 "gdT" = c('TRDC'), 
                                 "MAIT" = c('SLC4A10'),
                                 "NK T" = c('FGFBP2'),
                                 "NK" = c('KLRF1'), 
                                 "ILCs" = c('AREG') ,
                                 "GC B" = c('CD79A', 'BCL6'),
                                 "Naive B" = c('CD19', 'MS4A1', 'IGHD'),
                                 "Memory B" = c('CD27', 'CD83'),
                                 "IgA plasma" = c('IGHA1'),
                                 "IgG plasma" = c('IGHG1'), 
                                 "Monocyte" = c( 'C1QA', 'S100A8'),
                                 "Macrophage" = c('CD14', 'CD163', 'MMP12', 'C1QA'),
                                 "cDC1" = c('CLEC9A'),
                                 "cDC2" = c('CD1C'),
                                 "Lymphoid DC" = c('LAMP3'),
                                 "Mast" = c('CPA3'),
                                 "Neutrophil" = c('FCGR3B'),
                                 "Stromal-1" = c('IGFBP7', 'ADAMDEC1'), 
                                 "Stromal-2" = c('NRG1'),  
                                 "Stromal-3" = c('SOX6'), 
                                 "Myofibroblast" = c('ACTA2'), 
                                 "Arterial" = c('PLVAP'), 
                                 "Venous" = c('ACKR1'), 
                                 "Pericyte" = c('NOTCH3'), 
                                 "Contractile pericyte" = c('RERGL'), 
                                 "Smooth muscle" = c('DES'), 
                                 "Lymphatic endothelium" = c('CCL21'), 
                                 "Glial" = c('S100B')
)


atac.unionHeatmap.markers <- list("Stem" = c('EPCAM','LGR5'), 
                                 "EC1-1" = c('GSTA1'),
                                 "EC1-2" = c('APOA4'),
                                 "EC2-1" = c('CD24'),
                                 "EC2-2" = c('CEACAM5'),
                                 "Goblet1" = c('MUC2'),
                                 "Goblet2" = c('BEST2'),
                                 "BEST4" = c('BEST4'),
                                 "Paneth" = c('DEFA5'),
                                 "EEC" = c('CHGA'),
                                 "Tuft" = c('FYB1'),
                                 "CD4 Tcm" = c('CD3D', 'TCF7'), 
                                 "Treg" = c('CTLA4', 'FOXP3'), 
                                 "CD103- CD4 Trm" = c('IL7R'), 
                                 "CD103+ CD4 Trm" = c('ITGAE'), 
                                 "CD103+ CD8 Trm" = c('CD8A'), 
                                 "KLRG1+ CD8 Trm" = c('KLRG1'), 
                                 "gdT" = c('TRDC'), 
                                 "MAIT" = c('SLC4A10'),
                                 "NK" = c('KLRF1'), 
                                 "ILCs" = c('AREG') ,
                                 "GC B" = c('CD79A', 'BCL6'),
                                 "Naive B" = c('CD19', 'MS4A1', 'IGHD'),
                                 "Memory B" = c('CD27', 'CD83'),
                                 "Plasma cell" = c('XBP1'),
                                 "Monocyte" = c( 'C1QA', 'S100A8'),
                                 "Macrophage" = c('CD14', 'CD163', 'MMP12', 'C1QA'),
                                 "DC" = c('CLEC9A', 'CD1C'),
                                 "Mast" = c('CPA3'),
                                 "Neutrophil" = c('FCGR3B'),
                                 "Stromal-1" = c('IGFBP7', 'ADAMDEC1'), 
                                 "Stromal-2" = c('NRG1'),  
                                 "Stromal-3" = c('SOX6'), 
                                 "Myofibroblast" = c('ACTA2'), 
                                 "Arterial" = c('PLVAP'), 
                                 "Venous" = c('ACKR1'), 
                                 "Pericyte" = c('NOTCH3'), 
                                 "Lymphatic endothelium" = c('CCL21'), 
                                 "Glial" = c('S100B')
)






















