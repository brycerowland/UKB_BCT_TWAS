#!/bin/bash

out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/regenie_step1
cov_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/new_covs
scr_out=/pine/scr/b/r/bryce38/twas/UKB_temp/data/regenie_temp_out
log_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/out_files

while read chr gene_name pheno; do

sbatch -J UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step1 --mem 12GB -t 7-00:00:00 -n 8 --mail-type END --mail-user bryce.rowland@unc.edu --out ${log_dir}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step1.out --wrap "/proj/yunligrp/share/TWAS/regenie/regenie \
  --step 1 \
  --bed /proj/yunligrp/users/bryce/twas/UKB/data/minority_TWAS_replication/UKB_LD_pruned_genotypes/ukb_allChrs_LD_pruned_MAF0.01 \
  --covarFile ${cov_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_regenie_covFile \
  --phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
  --phenoCol $pheno \
  --bsize 1000 \
  --lowmem \
  --loocv \
  --lowmem-prefix ${scr_out}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno} \
  --out ${out_dir}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step1"

# /proj/yunligrp/share/TWAS/regenie/regenie \
#   --step 1 \
#   --bed /proj/yunligrp/users/bryce/twas/UKB/data/minority_TWAS_replication/UKB_LD_pruned_genotypes/ukb_allChrs_LD_pruned_MAF0.01 \
#   --covarFile ${cov_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_regenie_covFile \
#   --phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
#   --phenoCol $pheno \
#   --bsize 1000 \
#   --lowmem \
#   --loocv \
#   --lowmem-prefix ${scr_out}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno} \
#   --out ${out_dir}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step1

done < <( awk 'NR==FNR { a[$1"_"$2]; next} $2"_"$3 in a' 03_rerun_list <( cut -f 1,4,5 novel_by_class_genes_b37_pos.tsv ) )
