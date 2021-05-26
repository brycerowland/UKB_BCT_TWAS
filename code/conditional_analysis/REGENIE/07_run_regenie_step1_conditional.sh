#!/bin/bash

out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/analysis_results
cov_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/regenie_covariates
scr_out=/pine/scr/b/r/bryce38/twas/UKB_temp/data/regenie_temp_out
: '
for chr in 1; do
	for pheno_class in PLT; do
		#Many jobs failed after 3 days when running on all phenotypes. Going to try 5 days, but reduce the phenotypes considered in each analysis. 

		if [ "$pheno_class" -eq "WBC" ]; then
			pheno_list="white_blood_cell_count,lymphocyte_count,monocyte_count,neutrophil_count,eosinophil_count,basophil_count,lymphocyte_percentage,monocyte_percentage,neutrophil_percentage,eosinophil_percentage,basophil_percentage"
		elif [ "$pheno_class" -eq "PLT" ]; then
			pheno_list="platelet_count,mean_platelet_volume,plateletcrit,platelet_distribution_width"
		else
			pheno_list="red_blood_cell_count,hemoglobin_concentration,hematocrit_percentage,mean_cell_volume,mean_cell_hemoglobin,mean_cell_hemoglobin_concentration,red_blood_cell_width,reticulocyte_percentage,reticulocyte_count,mean_reticulocyte_volume,mean_sphered_cell_volume,immature_reticulocyte_fraction,hls_reticulocyte_percentage,hls_reticulocyte_count"
		fi

		sbatch -J UKB_TWAS_conditional_test_chr${chr}_${pheno_class}_regenie_step1 --mem 12GB -t 3-00:00:00 -n 8 --mail-type END --mail-user bryce.rowland@unc.edu --out ${out_dir}/UKB_TWAS_conditional_test_chr${chr}_${pheno_class}_regenie_step1.out --wrap "/proj/yunligrp/share/TWAS/regenie/regenie \
  --step 1 \
  --bed /proj/yunligrp/users/bryce/twas/UKB/data/minority_TWAS_replication/UKB_LD_pruned_genotypes/ukb_allChrs_LD_pruned_MAF0.01 \
  --covarFile ${cov_dir}/UKB_TWAS_conditional_test_covFile_chr${chr}_${pheno_class} \
  --phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
  --phenoColList $pheno_list \
  --bsize 1000 \
  --lowmem \
  --loocv \
  --lowmem-prefix ${scr_out}/UKB_TWAS_conditional_test_chr${chr}_${pheno_class} \
  --out ${out_dir}/UKB_TWAS_conditional_test_chr${chr}_${pheno_class}_regenie_step1"
	
	done
done
'

while read chr pheno_class; do
	if [ "$pheno_class" = "WBC" ]; then
                        pheno_list="white_blood_cell_count,lymphocyte_count,monocyte_count,neutrophil_count,eosinophil_count,basophil_count,lymphocyte_percentage,monocyte_percentage,neutrophil_percentage,eosinophil_percentage,basophil_percentage"
                elif [ "$pheno_class" = "PLT" ]; then
                        pheno_list="platelet_count,mean_platelet_volume,plateletcrit,platelet_distribution_width"
                else
                        pheno_list="red_blood_cell_count,hemoglobin_concentration,hematocrit_percentage,mean_cell_volume,mean_cell_hemoglobin,mean_cell_hemoglobin_concentration,red_blood_cell_width,reticulocyte_percentage,reticulocyte_count,mean_reticulocyte_volume,mean_sphered_cell_volume,immature_reticulocyte_fraction,hls_reticulocyte_percentage,hls_reticulocyte_count"
                fi
                sbatch -J UKB_TWAS_conditional_test_chr${chr}_${pheno_class}_regenie_step1 --mem 12GB -t 7-00:00:00 -n 8 --mail-type END --mail-user bryce.rowland@unc.edu --out ${out_dir}/UKB_TWAS_conditional_test_chr${chr}_${pheno_class}_regenie_step1.out --wrap "/proj/yunligrp/share/TWAS/regenie/regenie \
  --step 1 \
  --bed /proj/yunligrp/users/bryce/twas/UKB/data/minority_TWAS_replication/UKB_LD_pruned_genotypes/ukb_allChrs_LD_pruned_MAF0.01 \
  --covarFile ${cov_dir}/UKB_TWAS_conditional_test_covFile_chr${chr}_${pheno_class} \
  --phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
  --phenoColList $pheno_list \
  --bsize 1000 \
  --lowmem \
  --loocv \
  --lowmem-prefix ${scr_out}/UKB_TWAS_conditional_test_chr${chr}_${pheno_class} \
  --out ${out_dir}/UKB_TWAS_conditional_test_chr${chr}_${pheno_class}_regenie_step1"
done < 07_rerun_list

