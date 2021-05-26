#!/bin/bash

dosage_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/UKB_QC_LDPruned0.8_dosage_files
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/regenie_covariates

pheno_class=$1
chr=$2

awk '{ $2="snp"(NR-1); print $0 }' <( zcat ${dosage_dir}/UKB_QC_LDPruned0.8_chr${chr}_${pheno_class}.dosage.gz ) | cut -f 2,7- -d ' ' | tawk | awk '{ print $1, $0 }' | sed 's/snp0 snp0/FID IID/g' > ${out_dir}/UKB_TWAS_conditional_test_covFile_chr${chr}_${pheno_class}

