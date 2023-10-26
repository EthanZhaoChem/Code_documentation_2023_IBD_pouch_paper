library(sys)
library(dplyr)
library(stringr)



# 1. write and read sample IDs in this project
SampleIDs <- c("HA01-AC", "HA01-TI", "HA02-AC", "HA02-TI", "HA04-AC", "HA04-TI", "HA50AC", "HA50TI", "HA51AC", "HA51TI", "HA55AC", "HA55TI",
               "OR101-POU", "OR101-PP", "OR102-POU", "OR102-PP", "OR109-POU", "OR109-PP", "OR43-POU", "OR43-PP", "OR48-POU", "OR48-PP", "OR72POU", "OR72PP")
write.table(SampleIDs, '~/yuzhao1/work/final_RC2rna/metadata/SampleIDs.csv',
            col.names = 'SampleID', row.names = F)

SampleIDs_pppou <- c("OR101-POU", "OR101-PP", "OR102-POU", "OR102-PP", "OR109-POU", "OR109-PP", "OR43-POU", "OR43-PP", "OR48-POU", "OR48-PP", 
               "OR72POU", "OR72PP", 'OR53-POU', 'OR53-PP')
write.table(SampleIDs_pppou, '~/yuzhao1/work/final_RC2rna/metadata/SampleIDs_pppou.csv',
            col.names = 'SampleID', row.names = F)


SampleIDs <- read.csv("~/yuzhao1/work/final_RC2rna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID



# 2. cellranger
fastq.folders <- list.dirs('/home/yuzhao1/gca/Pouch_scRNA/fastq/', recursive = F, full.names = F)
problematic.samples <- c()

# cellranger commands
for(i in 1:length(fastq.folders)){
  i.fastq.folder <- fastq.folders[[i]]
  i.fastq.folder.path <- paste0('/home/yuzhao1/gca/Pouch_scRNA/fastq/', i.fastq.folder)
  i.fastq.folder.files <- list.files(path = i.fastq.folder.path, recursive = F)
  i.fastq.folder.samples <- c()
  
  # get unique sample names from a patient folder
  for (j.file in i.fastq.folder.files) {
    lane.start <- j.file %>% gregexpr('L00', .) %>% unlist()
    
    # identify samples that can't be grouped by lanes
    if(lane.start < 1){
      problematic.samples <- unique(c(problematic.samples, j.file))
      next
    }
    
    # identify the existed sample name
    sample.name <- substr(j.file, 1, lane.start-2)
    i.fastq.folder.samples <- unique(c(i.fastq.folder.samples, sample.name))
  }
  
  # print commands to submit job for each sample
  for (k.sample in i.fastq.folder.samples) {
    # clean up sampleID
    SampleID.start <- max(k.sample %>% gregexpr('HA', .) %>% unlist(),
                          k.sample %>% gregexpr('OR', .) %>% unlist())
    temp1 <- substr(k.sample, SampleID.start, nchar(k.sample))
    SampleID <- strsplit(temp1, split = '_S')[[1]][[1]]
    
    # remove S string in k.sample
    samplename.end <- as.numeric(str_locate(k.sample, SampleID))[2]
    k.sample.updated <- substr(k.sample, 1, samplename.end)
    cat('sbatch -J ', SampleID, ' cellranger.batch ',
        SampleID, " /home/yuzhao1/gca/Pouch_scRNA/fastq/", i.fastq.folder, ' ', k.sample.updated,
        '\n',sep = '')
  }
}

# print samples that can't be grouped by lanes
for (i.problem in problematic.samples) {
  cat(i.problem, '\n')
}

# check cellranger out slurm file, successful or not
SampleIDs <- read.csv("~/yuzhao1/work/final_RC2rna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID
successes <- 0
for (i in 1:length(SampleIDs)){
  
  # no need to re-process AC-TI data
  if(grepl('TI', SampleIDs[[i]]) | grepl('AC', SampleIDs[[i]])){
    next
  }
  i.outfile <- paste0('~/yuzhao1/work/final_RC2rna/upstream/', SampleIDs[[i]], '.out')
  temp <- readLines(i.outfile)
  flag <- grepl('Pipestance completed successfully!', temp) %>% sum()
  if(flag >0){
    successes <- successes+1
  }
  if(flag < 1){
    cat(i.outfile, '\n')
  }
}

# check the number of successfull samples is same with #samples: all successes, mv output logs to log folder
library(filesstrings)
cat(successes)
# 12 healthy samples and s inf samples(manually checked) are successfully finished
if(successes == 12){
  SampleIDs <- read.csv("yuzhao1/work/final_RC2rna/metadata/SampleIDs.csv")
  SampleIDs <- SampleIDs$SampleID
  
  SampleIDs <- SampleIDs[c(grep('POU', SampleIDs), grep('PP', SampleIDs))]
  out.files <- paste0('~/yuzhao1/work/final_RC2rna/upstream/', SampleIDs, '.out')
  err.files <- paste0('~/yuzhao1/work/final_RC2rna/upstream/', SampleIDs, '.err')
  for (log in c(out.files, err.files)){
    file.move(log, '~/yuzhao1/work/final_RC2rna/upstream/log_cellranger')
  }
}
length(list.files('~/yuzhao1/work/final_RC2rna/upstream/log_cellranger', recursive = F))



# 3. velocyto commands
SampleIDs <- read.csv("~/yuzhao1/work/final_RC2rna/metadata/SampleIDs_pppou.csv")
SampleIDs <- SampleIDs$SampleID
for (SampleID in SampleIDs){
  # no need to re-process AC-TI data
  if(grepl('TI', SampleID) | grepl('AC', SampleID)){
    next
  }
  
  sample.dir <- paste0('~/yuzhao1/work/final_RC2rna/upstream/', SampleID)
  
  cat('sbatch -J ', SampleID,
  ' velocyto.batch ',sample.dir,'\n')
}

# check velocyto out slurm file, successful or not
out.files <- paste0('~/yuzhao1/work/final_RC2rna/upstream/', SampleIDs, '.out')
successes <- 0
for (i.outfile in out.files) {
  temp <- readLines(i.outfile)
  flag <- grepl('Terminated Succesfully!', temp) %>% sum()
  if(flag >0){
    successes <- successes+1
  }
  if(flag < 1){
    cat(i.outfile, '\n')
  }
}

# check the number of successfull samples is same with #samples: all successes, mv output logs to log folder
library(filesstrings)
cat(successes)
if(successes == length(SampleIDs)){
  out.files <- paste0('~/yuzhao1/work/final_RC2rna/upstream/', SampleIDs, '.out')
  err.files <- paste0('~/yuzhao1/work/final_RC2rna/upstream/', SampleIDs, '.err')
  for (log in c(out.files, err.files)){
    file.move(log, '~/yuzhao1/work/final_RC2rna/upstream/log_velo/')
  }
}
length(list.files('~/yuzhao1/work/final_RC2rna/upstream/log_velo/', recursive = F))
















