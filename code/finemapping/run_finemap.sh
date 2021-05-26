#!/bin/bash

pheno=$1
locus_name=$2

finemap_dir=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/FINEMAP/
twas_results=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results

mkdir -p ${finemap_dir}/input_files/${pheno}
mkdir -p ${finemap_dir}/master_files/${pheno}
mkdir -p ${finemap_dir}/output_files/${pheno}

input_dir=${finemap_dir}/input_files/${pheno}
master_dir=${finemap_dir}/master_files/${pheno}
out_dir=${finemap_dir}/output_files/${pheno}
cor_dir=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/predicted_expression_correlation_matrices/${pheno}

# This file creates the necessary input to run FINEMAP on a given TWAS locus. 

# 1. Create the "z" file which contains TWAS summary statistics.

#Get list of genes at the locus
awk -v locus_name=$locus_name '$10==locus_name { print $4 }' ${twas_results}/UKB_BCT_TWAS_results_${pheno}_ref_file_TWASLoci > ${input_dir}/finemap_input_${locus_name}_gene_list

#Then pull their summary statistics and format
awk 'NR==FNR {a[$1]; next} $3 in a { print $3, $1, $2, "A", "T", 0.1, $8, $9}' ${input_dir}/finemap_input_${locus_name}_gene_list ${twas_results}/UKB_BCT_TWAS_results_${pheno}.regenie  | sed '1 i\rsid chromosome position allele1 allele2 maf beta se' > ${input_dir}/finemap_input_${locus_name}.z

#Remove header from corMatrix
tail -n +2 ${cor_dir}/${pheno}_${locus_name}_PredExp_corMatrix.ld > ${cor_dir}/${pheno}_${locus_name}_PredExp_corMatrix_FineMap.ld

rm -f ${master_dir}/finemap_master_${locus_name}

#Create finemap master file
echo "z;ld;snp;config;cred;log;n_samples" >> ${master_dir}/finemap_master_${locus_name}
echo "${input_dir}/finemap_input_${locus_name}.z;${cor_dir}/${pheno}_${locus_name}_PredExp_corMatrix_FineMap.ld;${out_dir}/${pheno}_${locus_name}.snp;${out_dir}/${pheno}_${locus_name}.config;${out_dir}/${pheno}_${locus_name}.cred;${out_dir}/${pheno}_${locus_name}.log;399835" >> ${master_dir}/finemap_master_${locus_name}

#Check if the number of genes at the locus is less than 5.
max_causal_genes=5
n_genes=$( cat ${input_dir}/finemap_input_${locus_name}_gene_list | wc -l )
if [ $n_genes -lt 5 ]; then
	max_causal_genes=$n_genes	
fi
#Run Finemap
../finemap_v1.4_x86_64/finemap_v1.4_x86_64 --sss --in-files ${master_dir}/finemap_master_${locus_name} --dataset 1 --n-causal-snps $max_causal_genes
