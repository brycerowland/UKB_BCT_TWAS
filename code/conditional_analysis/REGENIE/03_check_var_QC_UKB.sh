#!/bin/bash

snp_stat_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/UKB_snpStats_reported_varlists
GWAS_file=/proj/yunligrp/users/bryce/twas/UKB/data/BCX_GWAS_2020_known_variants/Vuckovic_SuppTables_clean.tsv
reported_var_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/UKB_QC_reported_vars

for chr in {22..1}; do
	for pheno_class in WBC PLT RBC; do
		echo $chr $pheno_class
		awk -v pheno_class=$pheno_class 'NR==FNR {a[$1]; next} $1=="associated_blood_index" || ($3 in a && $2==pheno_class)' OFS="\t" <( zcat ${snp_stat_dir}/UKB_SnpStats_chr${chr}_BCX_variants.txt.gz | grep -v ^# | awk -v chr=$chr '$1~alternate || ($10 > 20 && $11 > 20  && $18 > .4) { $1=chr":"$4"_"$5"_"$6; print $0 }' ) $GWAS_file > ${reported_var_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}
	done
done
