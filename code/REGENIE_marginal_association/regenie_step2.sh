#!/bin/bash

results_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/out_files

#Running step 2 of regenie, 

#Checking again with their paper, they used 64GB of memory, and 16 CPU's, but the analysis itself only took 24GB of memory.
sbatch -J UKB_BCT_REGENIE_step2 --mem 4GB -t 4:00:00 -n 8 --mail-type END --mail-user bryce.rowland@unc.edu --out ${out_dir}/UKB_BCT_REGENIE_step2.out \
--wrap "/proj/yunligrp/share/TWAS/regenie/regenie \
  --step 2 \
  --pgen /proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/predicted_gexp/DGN_UKB_intersection_allChr_trickDosage \
  --phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
  --bsize 1000 \
  --pred /proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/out_files/UKB_BCT_step1_pred.list \
  --split \
  --loocv \
  --out ${results_dir}/UKB_BCT_TWAS_results"
