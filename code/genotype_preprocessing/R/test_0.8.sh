#!/bin/bash

module add r/3.4.1

script_chr=$1

scr_data="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data"
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/model_training
UKB_snp_dir=/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data/UKB_b38_liftover_output

while read chr start_pos end_pos gene_name chunk
do

	#Assemble the snplist
	sed '1d' ${scr_data}/DGN_UKB_filtered_imputed_data/DC_output/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.mach.info | cut -f1 | awk 'BEGIN {FS=":"; OFS="\t"} {print chr$1, $2, $3, $4}' > ${scr_data}/DGN_UKB_filtered_imputed_data/snp_lists/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection_snpList 

	if [ ! -f ${UKB_snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_UKBOrientation ]; then
		cat <( awk -F ":" '{ print $1, $2, $4, $3 }' OFS=":" ${UKB_snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_Flipped ) ${UKB_snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_nonFlipped > ${UKB_snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_UKBOrientation
	fi

        Rscript EN_0.8.setseed.R \
		-d ${scr_data}/DGN_UKB_filtered_imputed_data/training_dosage/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_PredictionFormat.gz \
		-e /pine/scr/m/x/munan/DGN/expression_cis/exp.${gene_name}.txt \
		-s ${scr_data}/DGN_UKB_filtered_imputed_data/snp_lists/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection_snpList \
		-w 1e6 -b $start_pos -f $end_pos \
		-o output_${gene_name} \
		-r 13 \
		--var_list ${UKB_snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_UKBOrientation


done < <( head -10 /proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/gene_ref_files/chr${script_chr}/gencode_build38_chr${script_chr}_genes.txt   )
		#-o ${out_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection \
