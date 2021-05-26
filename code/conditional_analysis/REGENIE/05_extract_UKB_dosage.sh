#!/bin/bash

module add qctool

pheno_class=$1
chr=$2

data_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE
reported_var_dir=${data_dir}/UKB_QC_LDPruned0.8_reported_vars
output_varFile_dir=${data_dir}/UKB_QC_LDPruned0.8_dosage_files


awk 'NR>1 { print $5":"$6":"$7":"$8, "NA", $5,$6,$7,$8 }' ${reported_var_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_LDPruned0.8 | sort -k1,1 -Vu | awk '($3<10){ $1=0$1; $3=0$3 } {print $0}' | sed 's/\t/ /g' | sed '1 i\SNPID rsid chromosome position alleleA alleleB' > ${output_varFile_dir}/UKB_QC_LDPruned0.8_${pheno_class}_chr${chr}_varFile.txt

qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample -compare-variants-by 'position,alleles' -incl-variants ${output_varFile_dir}/UKB_QC_LDPruned0.8_${pheno_class}_chr${chr}_varFile.txt -incl-samples ~/scr/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -og ${output_varFile_dir}/UKB_QC_LDPruned0.8_chr${chr}_${pheno_class}.dosage.gz



