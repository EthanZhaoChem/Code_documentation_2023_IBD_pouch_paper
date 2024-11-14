dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(Seurat)
library(magrittr)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/manu/rc2/scripts/umap_colors.R')

################ read metadata ##############
out.dir <- '~/yuzhao1/work/final_RC2rna/0revision/reviewer2/plots/'
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

############### analyze a single gene ###############
gene <- 'NRG1'
p <- ggplot(df.healthy, aes_string(x='biopsy_location', y=gene, fill='biopsy_location',
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
                     label.y = c(1.9),
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = F)+
  facet_grid(. ~ disease_status, scales = "free_x", space = "free_x")+
  labs(x = "Biopsy location", y = "Normalized expression", title = gene)
  
pdf(paste0(out.dir, 'bulk_', gene, '.pdf'), width = 2, height = 5)
print(p)
dev.off()


gene <- 'CEACAM5'
p <- ggplot(df.healthy, aes_string(x='biopsy_location', y=gene, fill='biopsy_location',
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
                     label.y = c(10),
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = F)+
  facet_grid(. ~ disease_status, scales = "free_x", space = "free_x")+
  labs(x = "Biopsy location", y = "Normalized expression", title = gene)

pdf(paste0(out.dir, 'bulk_', gene, '.pdf'), width = 2, height = 5)
print(p)
dev.off()


############### analyze gene gene correlation ###############
# NRG1 and CEACAM5 in pouch
df.healthy.pouch <- df.healthy[df.healthy$biopsy_location=='Pouch',]
xx <- cor.test(df.healthy.pouch$CEACAM5, df.healthy.pouch$NRG1, method = 'spearman')
pvalue <- xx$p.value %>% format(., digits=2, scientific = T)
rho <- xx$estimate %>% format(., digits=2)

p <- ggplot(df.healthy.pouch, aes(x=CEACAM5, y=NRG1)) + 
  geom_smooth(method="lm", col="black", formula = 'y ~ x', fill = '#a6cee3', size=0.5) + 
  geom_point(size=0.5, color='#0570b0') + 
  theme_pubr() + # Use ggpubr theme
  theme(axis.title = element_text(size = 8, face = 'bold'),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "bottom") +
  labs(title = paste0('rho=', rho, ', p-value=', pvalue), x = "CEACAM5", y = "NRG1")


pdf(paste0(out.dir, 'bulk_cor_CEACAM5_NRG1_healthy_pouch.pdf'), height = 3, width = 3)
print(p)
dev.off()









