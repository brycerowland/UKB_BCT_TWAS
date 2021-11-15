#!/bin/bash

var_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/GWAS_variants_other_classes_within_1MB
new_covs=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/new_covs


dosage_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/UKB_QC_LDPruned0.8_dosage_files
old_cov_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/regenie_covariates

while read chr gene_name pheno_class; do

	echo $gene_name

	#Extract dosage values for any additional SNPs which were included in the SnpPickOutput file. 

	awk 'NR==FNR { a[$1]; next } $6 in a { print $2, $6 }' ${var_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_SnpPickOutput ${var_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class | sort -k2,2 -u > ${var_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_newCovs


	while read class; do
		awk 'NR==FNR {a[$1]; next} $4 in a' <( awk -v class=$class '$1==class { print $2}' ${var_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_newCovs ) <( zcat ${dosage_dir}/UKB_QC_LDPruned0.8_chr${chr}_${class}.dosage.gz ) 

	done < <( cut -f 1 -d ' ' ${var_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_newCovs | sort -u ) | awk '{$2="new_class_snp"(NR); print $0}' | cut -f 2,7- -d ' ' | tawk > ${new_covs}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_regenie_newCovs


	#paste to existing covariate files. 
	paste -d ' ' ${old_cov_dir}/UKB_TWAS_conditional_test_covFile_chr${chr}_${pheno_class} ${new_covs}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_regenie_newCovs > ${new_covs}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_regenie_covFile

done < <( cut -f 1,4,6 novel_by_class_genes_b37_pos.tsv | tail -n +2 | sort -u  )


