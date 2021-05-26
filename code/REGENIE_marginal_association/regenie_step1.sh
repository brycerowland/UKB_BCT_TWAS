#!/bin/bash

out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/out_files

#Running step 1 of regenie, which takes a genome-wide set of genotyped SNPS and phenotypes as input. 

#Checking again with their paper, they used 64GB of memory, and 16 CPU's, but the analysis itself only took 24GB of memory.
sbatch -J UKB_BCT_REGENIE_step1 --mem 12GB -t 3-00:00:00 -n 8 --mail-type END --mail-user bryce.rowland@unc.edu --out ${out_dir}/UKB_BCT_REGENIE_step1.out \
--wrap "/proj/yunligrp/share/TWAS/regenie/regenie \
  --step 1 \
  --bed /proj/yunligrp/users/bryce/twas/UKB/data/minority_TWAS_replication/UKB_LD_pruned_genotypes/ukb_allChrs_LD_pruned_MAF0.01 \
  --phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
  --bsize 1000 \
  --lowmem \
  --loocv \
  --lowmem-prefix ${out_dir}/temp_wbc \
  --out ${out_dir}/UKB_BCT_step1"
