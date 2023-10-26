dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')

################ read metadata ##############
out.dir <- '~/yuzhao1/work/manu/rc2/plots/extra_bulk_deconvolution/expression_boxPlot/'
metadata <- read.table('/project/gca/yuzhao1/work/final_RC2rna/bulk/metadata/metadatacurated.csv', sep = ',', header = T)
metadata$Patient_ID <- gsub('patient: ', '', metadata$Patient_ID)
metadata$Original_diagnosis <- gsub('diagnosis: ', '', metadata$Original_diagnosis)
metadata$disease_status <- gsub('prognosis: ', '', metadata$disease_status)
metadata$biopsy_location <- gsub('tissue: ', '', metadata$biopsy_location)
metadata$time <- gsub('biopsytime: ', '', metadata$time)
metadata$age <- gsub('age: ', '', metadata$age)
metadata$sex <- gsub('Sex: ', '', metadata$sex)
metadata$ethnicity <- gsub('ethnicity: ', '', metadata$ethnicity)
rownames(metadata) <- metadata$Sample_ID

################ prepare counts ##############
raw_data <- read.table('~/yuzhao1/work/final_RC2rna/bulk/dataset/GSE81266_Pouch.combat.txt',
                       sep = '\t', header = T)
counts <- raw_data[, -(1:4)]
rownames(counts) <- raw_data$Symbol
counts <- na.omit(counts)
counts <- t(counts)

################ get a union dataframe with metadata and counts ################
df <- merge(metadata, counts, by="row.names", all=TRUE) 
df$Row.names <- NULL

############### filter noise ###############
# remove FAP and Ileum samples because they are not informative
# remove sample EC11b, because it is sequenced twice, all other samples in same cohort are single end seq
df <- subset(df, biopsy_location!='Ile')
df <- subset(df, disease_status != 'FAP')
df <- subset(df, Sample_ID!='EC11b')

df$Patient_time <- paste0(df$Patient_ID, '-', df$time)
df$biopsy_location <- unlist(mapvalues(df$biopsy_location, 
                                           from = c("Ileal pouch", "Prepouch ileum"), 
                                           to = c('Pouch', 'PrePouch')))
df$time <- unlist(mapvalues(df$time, 
                                from = c("4mo", "8mo", "12mo", ">12mo"), 
                                to = c("4month", "8month", "12month", ">12month")))

############### analyze the df you need ###############
df.healthy <- subset(df, disease_status == 'Healthy') 
df.pouchitis <- subset(df, disease_status == 'Pouchitis') 

############### analyze a single gene: healthy ###############
gene <- 'CEACAM5'
p <- ggplot(df.healthy, aes(x=biopsy_location, y=CEACAM5, fill=biopsy_location,
                            outlier.shape = NA, 
                            bxp.errorbar=T, bxp.errorbar.width = 0.2)) + 
  geom_boxplot()+ 
  geom_point(position=position_jitterdodge())+
  # geom_jitter(color="black", size=1.2, alpha=0.9, width = 0.1) +
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 10, angle=45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),        legend.position = "none",
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        plot.title = element_text(hjust = 0.5)
  )+
  scale_fill_manual(values=c('PrePouch' = "#b3de69", 'Pouch' = "#fdb462"))+
  stat_compare_means(label.x = 1.5, 
                     label.y = c(9),
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = F)+
  facet_grid(. ~ time, scales = "free_x", space = "free_x")+
  labs(x = "Biopsy location", y = "Normalized expression", title = gene)+
  
  
pdf(paste0(out.dir, 'healthy_', gene, '.pdf'), width = 5, height = 5)
print(p)
dev.off()

############### analyze a single gene: pouchitis ###############
gene <- 'CEACAM5'
p <- ggplot(df.pouchitis, aes(x=biopsy_location, y=CEACAM5, fill=biopsy_location,
                            outlier.shape = NA, 
                            bxp.errorbar=T, bxp.errorbar.width = 0.2)) + 
  geom_boxplot()+ 
  geom_point(position=position_jitterdodge())+
  # geom_jitter(color="black", size=1.2, alpha=0.9, width = 0.1) +
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 10, angle=45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),        legend.position = "none",
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        plot.title = element_text(hjust = 0.5)
  )+
  scale_fill_manual(values=c('PrePouch' = "#b3de69", 'Pouch' = "#fdb462"))+
  stat_compare_means(label.x = 1.5, 
                     label.y = c(9),
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = F)+
  facet_grid(. ~ time, scales = "free_x", space = "free_x")+
  labs(x = "Biopsy location", y = "Normalized expression", title = gene)+
  
  
  pdf(paste0(out.dir, 'pouchitis_', gene, '.pdf'), width = 5, height = 5)
print(p)
dev.off()

################ heatmap healthy ##############
# get mean expression for each location~time
genes.all <- colnames(counts)
df.healthy$biopsy_location_time <- paste0(df.healthy$biopsy_location, ' ', df.healthy$time)
df.healthy.mean <- df.healthy %>% 
  group_by(biopsy_location_time) %>% 
  dplyr::summarize(across(all_of(genes.all), mean), .groups = 'drop')%>%
  as.data.frame()
df.healthy.mean$biopsy_location_time <- factor(df.healthy.mean$biopsy_location_time, 
                                        levels = )
# select a few genes to show pouch identity
interesting.genes <- c('REG1A', 'DMBT1', 'SPINK1', 'MUC2', 'OLFM4',
                       'ADH1C', 'SATB2', 'CEACAM7', 'CA2', 'AQP8', 
                       'HES1', 'NXPE4', 'HMGCS2','RPSA', 'CFTR', 'RPL11',
                       'CEACAM5', 'SOX9', 'MECOM', 'TFCP2L1', 'CD24','FCGBP',
                       'LCN2', 'SLC7A7', 'GSTA1', 'ACE2', 'RBP2',
                       'APOA4', 'FABP6', 'HLA−DRA','ABCC2', 'ENPEP',
                       'EFNA1', 'SAE1', 'APOA1','APOB', 'SELENOP', 'TRPM6')
interesting.genes %<>% intersect(., colnames(counts))

# build the data frame and plot
df_plot <- df.healthy.mean[, c('biopsy_location_time', interesting.genes)]
rownames(df_plot) <- df_plot$biopsy_location_time
df_plot$biopsy_location_time <- NULL
df_plot <- as.matrix(df_plot)
df_plot %<>% scale(.)
df_plot <- df_plot[c("Pouch 4month", "Pouch 8month", "Pouch 12month", "Pouch >12month",
                     "PrePouch 4month", "PrePouch 8month", "PrePouch 12month", "PrePouch >12month"), ]


pdf(paste0(out.dir, '00heatmap_healthy.pdf'), width = 12, height = 6)
col_fun = colorRamp2(seq(-2, 2, 4/8), rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')))
p1 <- Heatmap(df_plot,   
              col = col_fun, 
              row_split = factor(c(rep('Pouch', 4), c(rep('Pre Pouch', 4))),
                                    levels = c('Pouch', 'Pre Pouch')),
              row_gap = unit(5, "mm"), column_gap = unit(0, "mm"), 
              border = T,
              border_gp = gpar(col = "black", lty = 1, lwd = 1),
              
              column_title = NULL,
              cluster_columns = F, 
              show_column_dend = F,
              
              row_names_gp = gpar(fontsize = 12),
              cluster_rows = F, 
              show_row_dend = F, 
              show_row_names = T,
              heatmap_legend_param = list(title = "Mean expression z score"),
              use_raster = F)
print(p1)
dev.off()



################ heatmap pouchitis ##############
# get mean expression for each location~time
genes.all <- colnames(counts)
df.pouchitis$biopsy_location_time <- paste0(df.pouchitis$biopsy_location, ' ', df.pouchitis$time)
df.pouchitis.mean <- df.pouchitis %>% 
  group_by(biopsy_location_time) %>% 
  dplyr::summarize(across(all_of(genes.all), mean), .groups = 'drop')%>%
  as.data.frame()
df.pouchitis.mean$biopsy_location_time <- factor(df.pouchitis.mean$biopsy_location_time, 
                                               levels = )
# select a few genes to show pouch identity
interesting.genes <- c('REG1A', 'DMBT1', 'SPINK1', 'MUC2', 'OLFM4',
                       'ADH1C', 'SATB2', 'CEACAM7', 'CA2', 'AQP8', 
                       'HES1', 'NXPE4', 'HMGCS2','RPSA', 'CFTR', 'RPL11',
                       'CEACAM5', 'SOX9', 'MECOM', 'TFCP2L1', 'CD24','FCGBP',
                       'LCN2', 'SLC7A7', 'GSTA1', 'ACE2', 'RBP2',
                       'APOA4', 'FABP6', 'HLA−DRA','ABCC2', 'ENPEP',
                       'EFNA1', 'SAE1', 'APOA1','APOB', 'SELENOP', 'TRPM6')
interesting.genes %<>% intersect(., colnames(counts))

# build the data frame and plot
df_plot <- df.pouchitis.mean[, c('biopsy_location_time', interesting.genes)]
rownames(df_plot) <- df_plot$biopsy_location_time
df_plot$biopsy_location_time <- NULL
df_plot <- as.matrix(df_plot)
df_plot %<>% scale(.)
df_plot <- df_plot[c("Pouch 4month", "Pouch 8month", "Pouch 12month",
                     "PrePouch 4month", "PrePouch 8month", "PrePouch 12month"), ]


pdf(paste0(out.dir, '00heatmap_pouchitis.pdf'), width = 12, height = 6)
col_fun = colorRamp2(seq(-2, 2, 4/8), rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')))
p1 <- Heatmap(df_plot,   
              col = col_fun, 
              row_split = factor(c(rep('Pouch', 3), c(rep('Pre Pouch', 3))),
                                 levels = c('Pouch', 'Pre Pouch')),
              row_gap = unit(5, "mm"), column_gap = unit(0, "mm"), 
              border = T,
              border_gp = gpar(col = "black", lty = 1, lwd = 1),
              
              column_title = NULL,
              cluster_columns = F, 
              show_column_dend = F,
              
              row_names_gp = gpar(fontsize = 12),
              cluster_rows = F, 
              show_row_dend = F, 
              show_row_names = T,
              heatmap_legend_param = list(title = "Mean expression z score"),
              use_raster = F)
print(p1)
dev.off()



################ heatmap healthy: core TF validation ##############
# get mean expression for each location~time
genes.all <- colnames(counts)
df.healthy$biopsy_location_time <- paste0(df.healthy$biopsy_location, ' ', df.healthy$time)
df.healthy.mean <- df.healthy %>% 
  group_by(biopsy_location_time) %>% 
  dplyr::summarize(across(all_of(genes.all), mean), .groups = 'drop')%>%
  as.data.frame()
df.healthy.mean$biopsy_location_time <- factor(df.healthy.mean$biopsy_location_time, 
                                               levels = )
# select a few genes to show pouch identity
interesting.genes <- c('CDX1', 'FOXP1', 'KLF5', 'EHF',  'NFIA', 'GATA6',
                       'BACH1', 'ESRRG', 'MAF', 'TBX3',
                       'HNF4G',  'NR1H4', 'PPARA')
interesting.genes %<>% intersect(., colnames(counts))

# build the data frame and plot
df_plot <- df.healthy.mean[, c('biopsy_location_time', interesting.genes)]
rownames(df_plot) <- df_plot$biopsy_location_time
df_plot$biopsy_location_time <- NULL
df_plot <- as.matrix(df_plot)
df_plot %<>% scale(.)
df_plot <- df_plot[c("Pouch 4month", "Pouch 8month", "Pouch 12month", "Pouch >12month",
                     "PrePouch 4month", "PrePouch 8month", "PrePouch 12month", "PrePouch >12month"), ]


pdf(paste0(out.dir, '00heatmap_healthy_TFs.pdf'), width = 12, height = 6)
col_fun = colorRamp2(seq(-2, 2, 4/8), rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')))
p1 <- Heatmap(df_plot,   
              col = col_fun, 
              row_split = factor(c(rep('Pouch', 4), c(rep('Pre Pouch', 4))),
                                 levels = c('Pouch', 'Pre Pouch')),
              row_gap = unit(5, "mm"), column_gap = unit(0, "mm"), 
              border = T,
              border_gp = gpar(col = "black", lty = 1, lwd = 1),
              
              column_title = NULL,
              cluster_columns = F, 
              show_column_dend = F,
              
              row_names_gp = gpar(fontsize = 12),
              cluster_rows = F, 
              show_row_dend = F, 
              show_row_names = T,
              heatmap_legend_param = list(title = "Mean expression z score"),
              use_raster = F)
print(p1)
dev.off()



################ heatmap pouchitis: core TF ##############
# get mean expression for each location~time
genes.all <- colnames(counts)
df.pouchitis$biopsy_location_time <- paste0(df.pouchitis$biopsy_location, ' ', df.pouchitis$time)
df.pouchitis.mean <- df.pouchitis %>% 
  group_by(biopsy_location_time) %>% 
  dplyr::summarize(across(all_of(genes.all), mean), .groups = 'drop')%>%
  as.data.frame()
df.pouchitis.mean$biopsy_location_time <- factor(df.pouchitis.mean$biopsy_location_time, 
                                                 levels = )
# select a few genes to show pouch identity
interesting.genes <- c('CDX1', 'FOXP1', 'KLF5', 'EHF',  'NFIA', 'GATA6',
                       'BACH1', 'ESRRG', 'MAF', 'TBX3',
                       'HNF4G',  'NR1H4', 'PPARA')
interesting.genes %<>% intersect(., colnames(counts))

# build the data frame and plot
df_plot <- df.pouchitis.mean[, c('biopsy_location_time', interesting.genes)]
rownames(df_plot) <- df_plot$biopsy_location_time
df_plot$biopsy_location_time <- NULL
df_plot <- as.matrix(df_plot)
df_plot %<>% scale(.)
df_plot <- df_plot[c("Pouch 4month", "Pouch 8month", "Pouch 12month",
                     "PrePouch 4month", "PrePouch 8month", "PrePouch 12month"), ]


pdf(paste0(out.dir, '00heatmap_pouchitis_TFs.pdf'), width = 12, height = 6)
col_fun = colorRamp2(seq(-2, 2, 4/8), rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')))
p1 <- Heatmap(df_plot,   
              col = col_fun, 
              row_split = factor(c(rep('Pouch', 3), c(rep('Pre Pouch', 3))),
                                 levels = c('Pouch', 'Pre Pouch')),
              row_gap = unit(5, "mm"), column_gap = unit(0, "mm"), 
              border = T,
              border_gp = gpar(col = "black", lty = 1, lwd = 1),
              
              column_title = NULL,
              cluster_columns = F, 
              show_column_dend = F,
              
              row_names_gp = gpar(fontsize = 12),
              cluster_rows = F, 
              show_row_dend = F, 
              show_row_names = T,
              heatmap_legend_param = list(title = "Mean expression z score"),
              use_raster = F)
print(p1)
dev.off()
































