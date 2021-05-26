#!/bin/bash

pheno=$1

#Runs finemap for all twas loci for a phenotype. 

locus_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results

#Make TWAS loci list
cut -f 10 ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_ref_file_TWASLoci | awk '$1!="NA"' | tail -n +2 | sort  | uniq -c | awk '{ print $2, $1 }' OFS="\t"  > ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_locusInfo

while read locus_name; do
	../run_finemap.sh $pheno $locus_name	
done < <( awk '$2 > 1' ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_locusInfo )

#rm ${locus_dir}/${pheno}_locus_list
