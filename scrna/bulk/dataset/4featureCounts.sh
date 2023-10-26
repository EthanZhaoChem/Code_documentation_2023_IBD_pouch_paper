#!/bin/bash
ls -1 /project/gca/yuzhao1/work/final_RC2rna/bulk/dataset/sorted_bam/*sorted.bam > bam_list.txt # then separate into two: bam_list_singleEnd.txt and bam_list_pairedEnd.txt

bam_files_singleEnd=$(cat bam_list_singleEnd.txt | tr '\n' ' ')
featureCounts -T 10 -s 0  -t exon  -g gene_id -a /project/gca/yuzhao1/software/cellranger/refdata/refdata-gex-GRCh38-2020-A/genes/genes.gtf --extraAttributes gene_name --primary -M --ignoreDup -Q 10 --donotsort -o counts_matrix_singleEnd.txt $bam_files_singleEnd

bam_files_pairedEnd=$(cat bam_list_pairedEnd.txt | tr '\n' ' ')
featureCounts -T 10 -s 0  -t exon  -g gene_id -a /project/gca/yuzhao1/software/cellranger/refdata/refdata-gex-GRCh38-2020-A/genes/genes.gtf --extraAttributes gene_name --primary -M --ignoreDup -Q 10 --donotsort -o counts_matrix_pairedEnd.txt -B -P -p  --countReadPairs $bam_files_pairedEnd



