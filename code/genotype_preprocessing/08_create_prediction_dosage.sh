#!/bin/bash

chr=$1


data_dir=/proj/yunligrp/users/bryce/twas/UKB/data
dosage_dir=${data_dir}/prediction_files/subset_dosage/chr${chr}
gene_ref_file=${data_dir}/prediction_files/gene_ref_files/gencode_build38_allChr_TWAS_genes.txt
training_dir=/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data/model_training/chr${chr}
varID_dir=/proj/yunligrp/users/bryce/twas/UKB/data/model_training/chr${chr}/varIDs


for chunk in {1..42}; do

if  [ -f ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs.dosage.gz ]; then 
	#Make b37 varID's and liftOver to b38   
	if [ "${chr}" -lt 10 ]; then
		zcat ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs.dosage.gz | tail -n +2 | cut -f 1,4-6 -d ' ' | sed 's/ /:/g' | sed 's/^0//g' > ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs_ID_b37
	else
		zcat ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs.dosage.gz | tail -n +2 | cut -f 1,4-6 -d ' ' | sed 's/ /:/g' > ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs_ID_b37
	fi      

	#Use ID liftover
	~/bin/id_liftOver.sh ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs_ID_b37 ~/bin/chain_files/hg19ToHg38.over.chain ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs_ID_b38

	#Create prediction format file, with all variants. 

	echo Creating prediction file . . . 
	paste -d ' ' <( cat <( echo "snpID" )  ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs_ID_b38 ) <( zcat ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs.dosage.gz | cut -f 7- -d ' ' )  | gzip -c > ${dosage_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_modelSNPs_PredictionReady.dosage.gz

fi      
done

