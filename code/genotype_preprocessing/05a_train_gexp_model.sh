#!/bin/bash

module add r/3.6.0

chr=$1

scr_data="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data"
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/model_training
snp_dir=${scr_data}/DGN_UKB_filtered_imputed_data/snp_lists
DC_dir=${scr_data}/DGN_UKB_filtered_imputed_data/DC_output

while read start_pos end_pos gene_name chunk
do
	sed '1d' ${DC_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.mach.info | cut -f1 | awk 'BEGIN {FS=":"; OFS="\t"} {print chr$1, $2, $3, $4}' > ${snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_snpList 

	cat <( awk 'NR==FNR {a[$1]; next} $1 in a' <( tail -n +2 ${DC_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.mach.info | cut -f 1 ) ${scr_data}/UKB_b38_liftover_output/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_Flipped  | awk -F ":" '{ print $1, $2, $4, $3 }' OFS=":" ) \
		<( awk 'NR==FNR {a[$1]; next} $1 in a' <( tail -n +2 ${DC_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.mach.info | cut -f 1 ) ${scr_data}/UKB_b38_liftover_output/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_nonFlipped ) > ${snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_UKBOrientation

        Rscript /proj/yunligrp/users/bryce/twas/UKB/code/genotype_preprocessing/R/EN_0.8.setseed.R \
		-d ${scr_data}/DGN_UKB_filtered_imputed_data/training_dosage/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_PredictionFormat.gz \
		-e /pine/scr/m/x/munan/DGN/expression_cis/exp.${gene_name}.txt \
		-s ${snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_snpList \
		-w 1e6 -b $start_pos -f $end_pos \
		-o ${out_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection \
		-r 13 \
		--var_list ${snp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_UKBOrientation

done < <( cut -f 2-5 /proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/gene_ref_files/chr${chr}/gencode_build38_chr${chr}_genes.txt  )
