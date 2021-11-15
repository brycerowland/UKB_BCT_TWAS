#!/bin/bash

cov_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/new_covs
step1_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/regenie_step1
gexp_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/sig_gene_predicted_gexp
log_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/out_files
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/regenie_step2

while read chr gene_name pheno pheno_class; do

	sbatch -J UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step2 --mem 12GB -t 6:00:00 --mail-type END --mail-user bryce.rowland@unc.edu --out ${log_dir}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step2.out --wrap "/proj/yunligrp/share/TWAS/regenie/regenie \
		--step 2 \
		--pgen ${gexp_dir}/${pheno_class}_TWAS_bonferroniSigGenes_predGexp_chr${chr} \
		--phenoFile /proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv \
		--covarFile ${cov_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_regenie_covFile \
		--bsize 1000 \
		--pred ${step1_dir}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step1_pred.list \
		--phenoCol $pheno \
		--chr $chr \
		--out ${out_dir}/UKB_TWAS_conditional_test_novelByClass_${gene_name}_${pheno}_regenie_step2" 
		

done < <( cut -f 1,4,5,6 novel_by_class_genes_b37_pos.tsv  | tail -n +3  ) 
