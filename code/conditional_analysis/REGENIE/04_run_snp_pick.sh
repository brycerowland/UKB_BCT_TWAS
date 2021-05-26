#!/bin/bash

pheno_class=$1
chr=$2

data_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE
reported_var_dir=${data_dir}/UKB_QC_reported_vars
out_dir=${data_dir}/UKB_QC_LDPruned0.8_reported_vars

for chr in {1..22}; do
	for pheno_class in WBC PLT RBC; do	

	LD_file=../LD_files/LD_file_chr${chr}_filterRsq0.8.gz

	echo $chr $pheno_class . . .
#	cut -f 10 ${reported_var_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr} | tail -n +2 | sort | uniq | awk '$1!="NA"' > ${out_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_SnpPickInput
	
	#Run SnpPick
#	/proj/yunligrp/bin/snp-pick -s ${out_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_SnpPickInput -l $LD_file -t .8 -o ${out_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_SnpPickOutput

	#Output subset data
	awk 'NR==FNR {a[$1]; next} $10 in a || $1=="associated_blood_index"' ${out_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_SnpPickOutput ${reported_var_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}  > ${out_dir}/UKB_QC_BCX_reported_vars_${pheno_class}_chr${chr}_LDPruned0.8 

	done
done
