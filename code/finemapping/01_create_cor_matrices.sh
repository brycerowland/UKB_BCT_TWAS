#!/bin/bash

module add r/3.6.0
pheno=$1

locus_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results

#Make TWAS loci list
cut -f 10 ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_ref_file_TWASLoci | awk '$1!="NA"' | sort | uniq > ${locus_dir}/${pheno}_loci_list

mkdir -p /proj/yunligrp/users/bryce/twas/UKB/data/finemapping/predicted_expression_correlation_matrices/${pheno}
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/predicted_expression_correlation_matrices/${pheno}

while read locus_name; do		
	../R/01_create_cor_matrix.R --locus $locus_name --data ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_ref_file_TWASLoci -e ~/bryce_group/twas/UKB/data/prediction_files/predicted_gexp -o ${out_dir}/${pheno}_${locus_name}_PredExp_corMatrix
done <  ${locus_dir}/${pheno}_loci_list 

rm ${locus_dir}/${pheno}_loci_list

