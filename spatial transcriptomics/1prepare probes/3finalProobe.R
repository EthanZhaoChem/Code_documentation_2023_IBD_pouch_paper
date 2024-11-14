source('/project/gca/yuzhao1/work/manu/rc2/scripts/rna_annotation_markers.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')
source('~/yuzhao1/work/manu/rc2/scripts/tfs.R')

# 1. read basic information
probes_colon <- read.csv('~/yuzhao1/work/final_RC2rna/0revision/spatial/metadata/Xenium_hColon_preview_metadata.csv', header = T, skip = 1)
probes_colon <- probes_colon$Gene

labels1_deg <- c('DMBT1',  'MUC2', 'OLFM4',  'LCN2',
                 'ADH1C', 'SATB2', 'CEACAM7', 'CA2', 'AQP8', 'HES1', 'NXPE4', 'HMGCS2',
                 'RPSA', 'CFTR', 'RPL11', 'CEACAM5', 'SOX9', 'MECOM', 'TFCP2L1', 'CD24','FCGBP',
                 'GSTA1', 'ACE2', 'RBP2', 'APOA4', 'FABP6', 'HLA−DRA','ABCC2', 'ENPEP',
                 'SELENOP', 'TRPM6', 
                 "REG1B", "FOXD1", "GATA5", "KLF7", "REG4", "IRF8", 
                 "STAT6", "STAT5", "SATB1", 
                 "HOXB7","MT1A", "KLF11" ,"MT1A", "ABCB9",
                 "REG3A","TUBA1B",  "CELA3B",
                 "FABP1", "KRT20","ZBTB7B",
                 "DEFA6", "DEFA5", "CLCA1", "TFF1", 'ATOH1',
                 'GUCA2B')

labels2_ec_inf <- c(
  'REG4','DMBT1','SLC6A14','VNN1','TGM2','TFF1','HLA-DRB1','MGAT3','PI3','OLFM4','LPCAT1',
  'UNC5CL','C2','TRIM22','C4BPB','LCN2','MYEOV','CASP1','UNC13D','CFB','TRIM40','SESTD1','KDELR3',
  'LPIN1','ARFGAP3','TMEM92','S100A11', 'ELOVL7','GBP1','FFAR4','ITGA2','XDH','DUSP6','EFNA2','STAT1',
  'MAP3K5','SLC20A1','ZNF703','P2RX4','RNF186','MESP1','MUC4','SELENBP1','OASL','DEFB1','HMGCS2',
  'CA2','AQP8','ANPEP','CASD1','SPINK5','KAZN','SMAD9','HYI','GPX4','TNNC2')

labels3_goblet_inf <- c("ALDOB",
                        "REG4","XKR9","TIMP1","SERPINB5","TFF1","AQP3","S100P","LPCAT1","CEACAM6",
                        "C2CD4B","TRIM22","B3GALT5","LYZ","TFF2","SCD", "TGFBI", "KLK12", 
                        "SDR16C5", "PDIA4",  "MYEOV","LINC00261", "B3GNT7","SEC24D","GSTA1",  
                        "ANPEP",  "PRAP1","SPINK5","CASD1","RNASE1", 
                        "PRDX5","SATB2-AS1", "HMGCS2","CA2", "URAD","PPARG")

labels4_anno <- unlist(rna.unionHeatmap.markers)
names(labels4_anno) <- NULL

labels5_tf_genes_all <- c(
  "BACH1","BARX2","CDX1", "CDX2", "EHF","ESRRG","FOS","FOSL2","FOXO1",
  "FOXP1","GATA6","HBP1", "HNF4G","HOXB7","IRF9", "JDP2",
  "JUN","JUND", "KLF5", "MAF","MAFF", "MAZ",
  "MEF2D","NFE2L1", "NFIA", "NR1H4","NR1I2",
  "NR5A2","PITX2","PPARA","SOX4", "TBX3","VDR" 
)



genes.important <- unique(c(labels1_deg, labels2_ec_inf, labels3_goblet_inf, labels4_anno, labels5_tf_genes_all))

# 2. select genes based on colon panel
addon_genes <- sort(setdiff(genes.important, probes_colon))
genes_remove <- c(
  # unuseful inflamed genes in Smillie only
  'FOS', 'SLC6A14', 'TFF1', 'PI3', 'MYEOV', 'CASP1', 'UNC13D', 'LPIN1', 'FFAR4', 
  'EFNA2', 'STAT1', 'MAP3K5', 'KAZN', 'SMAD9', 'SLC20A1', 'P2RX4', 'TFF2', 'SCD', 
  'KLK12', 'SDR16C5', 'PDIA4', 'MYEOV', 'LINC00261', 'B3GNT7', 'SEC24D', 'ANPEP', 
  'URAD', 'PPARG',
  # unuseful TF
  'BARX2', 'FOXO1', 'HBP1', 'JDP2', 'MAZ', 'MEF2D', 'NFE2L1', 'NR1I2', 'NR5A2', 'VDR', 'GATA6',
  # unuseful cell type annotation
  'S100B', 'RERGL', 'ACKR1', 'PLVAP',  
  'LAMP3', 'CD83', 'MS4A1', 'CD79A', 'IGHD', 'IGHG1', 'IGHA1', 'TCF7', 'EPCAM', 'IGFBP7',
  # low priority genes
  'ZNF703', 'REG4', 'CFB', 'SESTD1', 'MAFF', 'S100A11', 'RPL11', 'OASL',
  'HLA−DRA', 'STAT5', 'HLA-DRB1', 'RPSA', 'DEFA5', 'MT1A', 'HOXB7', 'KRT20',
  'ABCB9', 'CELA3B', 'KLF7', 'KLF11', 'ZBTB7B', 'XDH', 'TUBA1B', 'XKR9', 'SATB2-AS1',
  'TMEM92', 'FABP1'
  )
addon_genes <- setdiff(addon_genes, genes_remove)

# 3. add important customized genes
labels7_customized <- c('GREM1', 'OGN', 'NRG1', 'NPY', 'PECAM1', 'SPI1', # annotation
                        'GATA4', 'PRDM1', 'FOSL1', 'FOSL2', 'IRF4', 'IRF9', 'HNF4A', # TF
                        'SATB2', 'SATB1', 'LRRC32', 'CIITA', 'CLDN2', 'CLDN15', 'PLAU') # interest
addon_genes <- sort(unique(c(addon_genes, labels7_customized)))

# 4. write to file
addon_genes
genes_alternative <- c("TGFBI", "TGM2", "TIMP1", "TNNC2", "TRDC", 'ACKR1', 'PLVAP',
                       "PITX2", "TRPM6", "UNC5CL", "VNN1", "SERPINB5")
writeLines(setdiff(addon_genes, genes_alternative), '~/yuzhao1/work/final_RC2rna/0revision/spatial/genes_add/genes_main.txt')
writeLines(genes_alternative, '~/yuzhao1/work/final_RC2rna/0revision/spatial/genes_add/genes_alternative.txt')







