#!/bin/bash

chr=$1

vcf_dir=/proj/yunligrp/users/jwen/copy_scr_munan/DGN/cohort_QC/DGN/chunking.raw/chunking.chr${chr}
DC_out_dir=/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data/DGN_unfiltered_dosages/chr${chr}


for chunk in {1..42}; do
	if [ -f ${vcf_dir}/DGN.chr${chr}.chunk${chunk}.vcf.gz ]; then
	#Run DC
        /proj/yunligrp/users/jdrosen/bin/DosageConvertor/bin/DosageConvertor --vcfDose ${vcf_dir}/DGN.chr${chr}.chunk${chunk}.vcf.gz \
        	--prefix ${DC_out_dir}/DGN_chr${chr}_chunk${chunk} \
        	--type mach \
        	--format 1
	fi
done
