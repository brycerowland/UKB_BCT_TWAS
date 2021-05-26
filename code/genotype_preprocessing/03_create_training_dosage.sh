#!/bin/bash

data_dir="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data/"
vcf_dir=${data_dir}/DGN_UKB_filtered_imputed_data/DGN_filtered_vcf
DC_out_dir=${data_dir}/DGN_UKB_filtered_imputed_data/DC_output
training_dosage_dir=${data_dir}/DGN_UKB_filtered_imputed_data/training_dosage

chr=$1

for chunk in {1..42}; do
	if [ -f ${vcf_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.vcf.gz ]; then
	#Run DC
        /proj/yunligrp/users/jdrosen/bin/DosageConvertor/bin/DosageConvertor --vcfDose ${vcf_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.vcf.gz \
        	--prefix ${DC_out_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection \
        	--type mach \
        	--format 1

        #Create training input dosage files
        paste <( zcat ${DC_out_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.mach.dose.gz | cut -f 1 | cut -f -1 -d '-' ) <( zcat ${DC_out_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.mach.dose.gz | cut -f 3- ) | gzip -c > ${training_dosage_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_PredictionFormat.gz
	fi
done
