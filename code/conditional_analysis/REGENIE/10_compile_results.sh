#!/bin/bash

results_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/analysis_results
marginal_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results 

for pheno in basophil_count basophil_percentage eosinophil_count eosinophil_percentage hematocrit_percentage hemoglobin_concentration hls_reticulocyte_count hls_reticulocyte_percentage immature_reticulocyte_fraction lymphocyte_count lymphocyte_percentage mean_cell_hemoglobin mean_cell_hemoglobin_concentration mean_cell_volume mean_platelet_volume mean_reticulocyte_volume mean_sphered_cell_volume monocyte_count monocyte_percentage neutrophil_count neutrophil_percentage platelet_count plateletcrit platelet_distribution_width red_blood_cell_count red_blood_cell_width reticulocyte_count reticulocyte_percentage white_blood_cell_count; do
	
	#combined results files
	cat <( head -1 ${results_dir}/byChr/UKB_BCT_TWAS_results_*_chr1_${pheno}.regenie | awk '{print $0, "PHENO"}' ) \
	<( for chr in {1..22}; do tail -n +2 ${results_dir}/byChr/UKB_BCT_TWAS_results_*_chr${chr}_${pheno}.regenie | awk -F " " -v pheno=$pheno '{$(NF+1)=pheno; print}'; done ) > ${results_dir}/UKB_BCT_TWAS_conditionalAnalysis_${pheno}_results

	#Ref file with gene info
	cat <( awk '{$(NF+1)="phenotype"; $(NF+1)="log10_conditional_p"; print }'  <( head -1 ${marginal_dir}/UKB_BCT_TWAS_results_${pheno}_ref_file_TWASLoci ) ) \
	<( awk -v pheno=$pheno 'NR==FNR {a[$3]=$11} $4 in a {$(NF+1)=pheno; $(NF+1)=a[$4]; print}' ${results_dir}/UKB_BCT_TWAS_conditionalAnalysis_${pheno}_results ${marginal_dir}/UKB_BCT_TWAS_results_${pheno}_ref_file_TWASLoci ) > ${results_dir}/UKB_BCT_TWAS_conditionalAnalysis_${pheno}_refFile

done 
