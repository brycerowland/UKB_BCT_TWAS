#!/bin/bash

var_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/UKB_QC_reported_vars
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/GWAS_variants_other_classes_within_1MB

old_snp_pick_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/UKB_QC_LDPruned0.8_reported_vars

#awk '{ $2="snp"(NR-1); print $0 }' <( zcat ${dosage_dir}/UKB_QC_LDPruned0.8_chr${chr}_${pheno_class}.dosage.gz ) | cut -f 2,7- -d ' ' | tawk | awk '{ print $1, $0 }' | sed 's/snp0 snp0/FID IID/g' > ${out_dir}/UKB_TWAS_conditional_test_covFile_chr${chr}_${pheno_class}

while read chr gene_name pheno_class start_pos end_pos; do

	LD_file=../../LD_files/LD_file_chr${chr}_filterRsq0.8.gz	

	echo $chr $gene_name $pheno_class 
	#first get variants within +/- 1MB which are not in the same phenotype class as the reported
	# novel TWAS association.
	awk -v class=$pheno_class -v start_pos=$start_pos -v end_pos=$end_pos '$2 != class && $1 !~ /assoc/ && $6 > (start_pos-1e6) && $6 < (end_pos+1e6)' ${var_dir}/UKB_QC_BCX_reported_vars_*_chr${chr} > ${out_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class


	#Create SNPPick Input
	cat ${old_snp_pick_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_SnpPickOutput <( cut -f 6 ${out_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class ) | sort -u > ${out_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_SnpPickInput

	#Run SNP Pick
	/proj/yunligrp/bin/snp-pick -s ${out_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_SnpPickInput -i ${old_snp_pick_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_SnpPickOutput -l $LD_file -t .8 -o ${out_dir}/chr${chr}_${gene_name}_nearby_GWAS_variants_other_BCT_class_SnpPickOutput


done < <( cut -f 1,4,6- novel_by_class_genes_b37_pos.tsv | tail -n +2 )

