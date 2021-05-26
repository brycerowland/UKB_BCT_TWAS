#!/bin/bash

data_dir=/proj/yunligrp/users/bryce/twas/UKB/data
results_dir=${data_dir}/conditional_analysis/REGENIE/analysis_results
gexp_dir=${data_dir}/conditional_analysis/REGENIE/sig_gene_predicted_gexp
cov_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/regenie_covariates
out_dir=${data_dir}/conditional_analysis/REGENIE/out_files

#Running step 2 of regenie, 

for chr in {22..1}; do
for pheno_class in RBC WBC PLT; do

 if [ "$pheno_class" = "WBC" ]; then
                        pheno_list="white_blood_cell_count,lymphocyte_count,monocyte_count,neutrophil_count,eosinophil_count,basophil_count,lymphocyte_percentage,monocyte_percentage,neutrophil_percentage,eosinophil_percentage,basophil_percentage"
                elif [ "$pheno_class" = "PLT" ]; then
                        pheno_list="platelet_count,mean_platelet_volume,plateletcrit,platelet_distribution_width"
                else
                        pheno_list="red_blood_cell_count,hemoglobin_concentration,hematocrit_percentage,mean_cell_volume,mean_cell_hemoglobin,mean_cell_hemoglobin_concentration,red_blood_cell_width,reticulocyte_percentage,reticulocyte_count,mean_reticulocyte_volume,mean_sphered_cell_volume,immature_reticulocyte_fraction,hls_reticulocyte_percentage,hls_reticulocyte_count"
                fi

#Checking again with their paper, they used 64GB of memory, and 16 CPU's, but the analysis itself only took 24GB of memory.
sbatch -J UKB_BCT_TWAS_conditionalAnalysis_chr${chr}_${pheno_class} -o ${out_dir}/UKB_BCT_TWAS_conditionalAnalysis_chr${chr}_${pheno_class}.out --mem 12GB -t 3:00:00 --mail-type END --mail-user bryce.rowland@unc.edu --wrap "/proj/yunligrp/share/TWAS/regenie/regenie \
  --step 2 \
  --pgen ${gexp_dir}/${pheno_class}_TWAS_bonferroniSigGenes_predGexp_chr${chr} \
  --phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
  --covarFile ${cov_dir}/UKB_TWAS_conditional_test_covFile_chr${chr}_${pheno_class} \
  --bsize 1000 \
  --pred ${results_dir}/UKB_TWAS_conditional_test_chr${chr}_${pheno_class}_regenie_step1_pred.list \
  --split \
  --phenoColList $pheno_list \
  --chr ${chr} \
  --out ${results_dir}/UKB_BCT_TWAS_results_${pheno_class}_chr${chr}"

done
done
