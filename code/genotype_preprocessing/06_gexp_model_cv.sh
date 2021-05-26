#!/bin/bash

module add r/3.4.1

script_chr=$1

scr_data="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data"
out_dir=${scr_data}/model_training

while read chr start_pos end_pos gene_name chunk
do
        sed '1d' ${scr_data}/DGN_UKB_filtered_imputed_data/DC_output/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.mach.info | cut -f1 | awk 'BEGIN {FS=":"; OFS="\t"} {print chr$1, $2, $3, $4}' > ${scr_data}/DGN_UKB_filtered_imputed_data/snp_lists/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection_snpList

        Rscript /proj/yunligrp/users/munan/scripts/EN_CV_0.1.setseed.R \
                -d ${scr_data}/DGN_UKB_filtered_imputed_data/training_dosage/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_PredictionFormat.gz \
                -e /proj/yunligrp/users/munan/DGN/imputed/LDsamples/expression_cis/exp.${gene_name}.txt \
                -s ${scr_data}/DGN_UKB_filtered_imputed_data/snp_lists/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection_snpList \
                -w 1e6 -b $start_pos -f $end_pos \
                -o ${out_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection \
                -r 13

done < /proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/gene_ref_files/chr${script_chr}/gencode_build38_chr${script_chr}_genes.txt
