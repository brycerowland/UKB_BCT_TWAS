#!/bin/bash

#This file performs gene expression imputation in the UKB European cohort for a given chromosome and 10Mb "chunk", numbered as DGN chunk in the gene_ref_file. 

chr=$1
chunk=$2

data_dir=/proj/yunligrp/users/bryce/twas/UKB/data
dosage_dir=${data_dir}/prediction_files/subset_dosage/chr${chr}
gene_ref_file=${data_dir}/prediction_files/gene_ref_files/gencode_build38_allChr_TWAS_genes.txt
training_dir=${data_dir}/model_training/chr${chr}
out_dir=${data_dir}/prediction_files/predicted_gexp/chr${chr}

module add r/3.6.0

while read start_pos end_pos gene_name; do
	echo $gene_name

        #Get variants from EN training file for the given gene
        cut -f 1-4 ${training_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection.betas.EN.txt  | tail -n +2 | sed 's/\t/:/g' | sed 's/^chr//g' > ${dosage_dir}/${gene_name}_EN_list

	#Pull dosages from UKB in order to do prediction. Modify as needed to fit prediction cohort. 
        awk 'NR==FNR {a[$1];next} $1 in a || $1=="snpID"' ${dosage_dir}/${gene_name}_EN_list <( zcat ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs_PredictionReady.dosage.gz ) | gzip -c  > ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection_modelSNPs_PredictionReady.dosage.gz

	#Run prediction
	Rscript ../pred_transpose_0.6.R \
		--dosage ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection_modelSNPs_PredictionReady.dosage.gz \
		-b ${training_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection.betas.EN.txt \
		-l ${training_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection.log \
		-o ${out_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection

	rm ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection_modelSNPs_PredictionReady.dosage.gz ${dosage_dir}/${gene_name}_EN_list

done < <( awk -v chr=$chr -v chunk=$chunk '$6==1 && $1==chr && $5==chunk' $gene_ref_file | cut -f 2-4  )
