#!/bin/bash

DGN_dir="/proj/yunligrp/users/jwen/copy_scr_munan/DGN/cohort_QC/DGN/chunking.raw"
out_path="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data/DGN_b38_var_lists"

module add samtools
chr=$1
vcf_dir=${DGN_dir}/chunking.chr${chr}

for chunk in {1..42}; do
	echo $chr $chunk ...
        lines=$( zcat ${vcf_dir}/DGN.chr${chr}.chunk${chunk}.vcf.gz | grep -v "^#" | wc -l )
        if [ "$lines" -ne "0" ]
        then
                # Filter for very well imputed variants | combine multiallelic variants | filter out multiallelic variants | get ID's
                bcftools filter -i "R2>0.8 && MAF>0.05" ${vcf_dir}/DGN.chr${chr}.chunk${chunk}.vcf.gz | bcftools norm -m +any | bcftools view -H --max-alleles 2 | cut -f 3 > ${out_path}/chr${chr}/DGN_b38_chr${chr}_chunk${chunk}_varIDs
        fi
done
