#!/bin/bash

GWAS_file=/proj/yunligrp/users/bryce/twas/UKB/data/BCX_GWAS_2020_known_variants/Vuckovic_SuppTables_clean.tsv

for pheno_class in WBC RBC PLT; do
	for chr in {1..22}; do
		out_file=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/reported_var_lists/${pheno_class}/BCX_chr${chr}_${pheno_class}_reported_vars
		awk -v chr=$chr -v pheno_class=$pheno_class 'NR==1 || (($5==chr) && ($2==pheno_class))' $GWAS_file > ${out_file}
	done
done
