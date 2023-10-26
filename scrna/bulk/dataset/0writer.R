
# fastq to sam, single end
for (i in 3493774:3493789){
  x <- paste0('sbatch -J SRR', i, ' 1fastq_to_sam_single_end.batch /project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/fastq/SRR', 
              i, '.fastq /project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/sam/SRR', i, '.sam\n')
  cat(x)
}

# fastq to sam, paired end
for (i in 3493790:3493850){
  x <- paste0('sbatch -J SRR', i, ' 1fastq_to_sam_paired_end.batch ',
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/fastq/SRR', i, '_1.fastq ',
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/fastq/SRR', i, '_2.fastq ',
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/sam/SRR', i, '.sam\n')
  cat(x)
}

# sam to sorted bam
for(i in 3493774:3493850){
  x <- paste0('sbatch -J ', i, ' 2sam_to_sortedBAM.batch ',
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/sam/SRR', i, '.sam ', 
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/bam/SRR', i, '.bam ', 
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/sorted_bam/SRR', i, '_sorted.bam\n'
  )
  cat(x)
}

# stringtie
for(i in 3493774:3493850){
  x <- paste0('sbatch -J ', i, ' 3stringtie.batch ',
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/sorted_bam/SRR', i, '_sorted.bam ',
              '/project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/gtf/SRR', i, '.gtf ', 
              'SRR', i, '\n'
  )
  cat(x)
}
















