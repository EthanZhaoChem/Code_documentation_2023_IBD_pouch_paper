################################# 1. ECs ###################################
## 1. ECs, five categories of genes:

# Category 1.1: TFs that affect the different ECs. (intersecting TF list and genes0.01, filter based on Vln plots)
rna_deg_markers_ec_tfs <- c('CREB3L3', 'MAF', 'SOX6', 'MYRFL', 'OSR2', 'TBX3', 'CREB3L2', 'CREB3L1', 'NR1H4', 'NR3C1', 'TFCP2L1',
                 'ZBTB10', 'KLF12', 'MYC', 'FOXD2', 'MYB', 'ZBTB7C', 'ARID5B', 'EHF', 'MECOM', 'HES1', 'SOX9', 'SATB2', 'PITX2',
                 'HMGA2', 'MLXIP')

# Category 1.2: ileum functions (from paper, TIvsAC, POU1vsPOU2, EC1vsEC2)
rna_deg_markers_ec_ileum <- c(
  'SLC10A2', 'FABP6','SLC51A','SLC51B', 'SLC5A1', 'SLC7A9', 'SLC15A1', 'SLC6A20',
  'SLC6A19', 'SLC7A7', 'SLC1A1', 'SLC3A1',
  
  'ABCC2', 'ABCG2', 'ABCB1',
  'APOA1', 'APOA4', 'APOB',
  'HLA-DRA', 'HLA-DRB1',

  'MTTP', 'MME', 'MGAM', 'MAF', 'MEP1B', 'MUC3A', 'MYO1A',
  'MOGAT2', 'PRAP1', 'VIL1', 'RBP2', 'DPEP1', 'ENPEP', 'SERPINA1', 'CREB3L3',
  'PMP22', 'XPNPEP2', 'GSTA1', 'DAB1', 'KHK', 'CYBRD1', 'ACE2', 'ACE', 'SOX6',
  'CIDEB', 'ALPI', 'TM4SF20', 'LIPA', 'DPP4', 'SELENOP', 'TM6SF2', 'ENTPD7',
  'FAM3C', 'SMIM24', 'DNASE1', 'GDA', 'CYP2C18', 'FBP1', 'GK', 'MS4A8', 'OAT',
  'PEPD', 'CDHR2', 'CLDN15', 'LPGAT1', 'DENND5B', 'SAT2', 'HEBP1', 'IL32', 'SCIN',
  'PCK2', 'TMIGD1', 'CA13', 'CNDP2', 'AIG1', 'TNFRSF1A', 'MEP1A', 'CBR1', 'BTNL8',
  
  # ileun only
  'TRPM6', 'G6PC', 'ENPP3', 'SMIM2-AS1', 'CHRNA7', 'GRAMD1C', 'PFKFB4', 'EFNA1',
  'SAE1', 'TRHDE', 'KCNJ3'
)

# Category 1.3: colon functions (from ACvsTI, POU2vsPOU1, EC2vsEC1)
rna_deg_markers_ec_colon <- c(
  #share
  'FCGBP', 'SLC38A1', 'MECOM', 'CD24', 'HES1', 'SLC12A2', 'RACK1', 'CFTR', 'PPIA',
  'GOLM1', 'PLCE1','EEF1G','PKP2','NACA','MGST1','MLEC','C1QBP','SNHG29','NME2',
  'CEACAM5','NUPR1','TFCP2L1','SLPI','PARM1','SLC39A8','NANS','SOX9','CXCL3',
  
  # RP genes
  "RPL22", "RPL11", "RPS8", "RPS7", "RPSA",

  #ac
  'AQP8', 'CEACAM7', 'MS4A12', 'CA2', 'SELENBP1', 'SATB2', 'HMGCS2', 'PRDX6', 'ITM2C',
  'MUC4', 'ADH1C', 'NXPE1', 'PCCA', 'GPX2', 'CA12', 'CKB', 'RBFOX1', 'PPP1R1B', 'NXPE4',
  'VSIG2', 'LEFTY1',  'IQCM', 'HMGA2',
  
  #pou2
  'DMBT1','REG1A', 'OLFM4','PITPNC1','ABCC1','MUC2', 'LCN2', 'SPINK1')
  
## 1. ECs, five categories of genes:
rna_deg_markers_stem_tfs <- c(
  "FOS","FOSB", "JUND", "SATB2","DACH1","ZBTB10","ZNF804A", "RORA", 
  "KLF4", "HES1", "PPARG","ZNF704","JUNB", "ELF3", "EHF","SOX6",
  "HNF4G","KLF3", "MECOM","GLIS3","TFCP2L1", "KLF12","MLXIPL","NR1H4",
  "NR5A2","NME2", "BNC2", "NFIB", "VDR",  "FOXP2","LCOR" 
)

rna_deg_markers_stem_nonTFs <- c(
"REG1B", "PLCG2", "LCN2", "MTRNR2L8", "MTRNR2L12", "REG1A", "SLC35F1", "AC105402.3", "UNC5D", "DMBT1", "ALDOB", 
"LAMA3", "RBFOX1", "MTRNR2L1", "SNHG5", "SPINK1", "MT-ND6", "GRIA4", "TEX14", "LYZ", "ROBO2", "CD74", 
"REG3A", "MIF", "TMEM238", "LHFPL3-AS2", "IER2", "PLA2G2A", "PLCB4", "KHDRBS2", "HSPA8", "BX284613.2", "HSPE1", 
"KIF26B", "OLFM4", "AGR2", "MT1G", "ADAMTSL1", "HMGCS2", "CCL25", "RPSA", "LGALS4", "PDE10A", "HSP90AA1", 
"CAMK2N1", "CADPS", "HSPH1", "CD9", "RHPN2", "DEFA5", "MALRD1", "EDIL3", "XACT", "C10orf99", "AC019330.1", 
"L1TD1", "CHRM3", "DEFA6", "SHISA9", "CPS1", "CPA6", "IFI27", "SI", "RAMP1", "IQCM", "KCNJ3", 
"TOX", "MAP3K20-AS1", "LDHB", "LEFTY1", "AL391117.1", "GALNT13", "PRDX5", "C15orf48", "ITM2C", "FUT8", "CNTN4", 
"DPP10", "XIST", "ENPEP", "LGALS2", "NUPR1", "CYP3A5", "MT1H", "HTR4", "CD24", "OAT", "RPS4Y1", 
"LRMDA", "FARP1", "UGT2B17", "SMOC2", "PRSS2", "IGFBP2", "MT2A", "HLA-B", "S100A6", "RIMS2", "CEACAM5", 
"PLXDC2", "FAM3D", "KRT20", "FREM1", "MUC4", "PDE3A"
)

rna_deg_markers_goblet <- c('FRY', 'MUC17', 'FAXDC2', #(direct target of FOXA3)
                            'SELENOP', #（direct target of SPDEF）
                            'INPP4B', 'GPRC5A', 'GPR39', 'ANO1', 'MGST1',
                            'CD24', 'MGST1', 'CTSC', 'IL18', 'CMSS1', 'NAP1L1',
                            'WFDC2', 'CCND1', 'DUT', 'SH3RF3', 'KCNJ3',
                            'IL1R2', 'MIR31HG', 'MGAM', 'MME', 'CCL25',
                            'CFTR', 'XPNPEP2', 'EMB', 'MPP6', 'SLC15A1', 'RIMS2', 
                            'CYBRD1', 'FREM1', 'ABCG2', 'DPEP1', 'RGS2', 'CYP4F3',
                            'GSTA1', 'L1TD1', 'MUC5B', 'LDHB')



####################################################################


