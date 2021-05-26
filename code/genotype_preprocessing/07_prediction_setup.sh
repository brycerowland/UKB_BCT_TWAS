#!/bin/bash

module add qctool

script_chr=$1
script_chunk=$2

bryce_data_dir="/proj/yunligrp/users/bryce/twas/UKB/data"
gencode_ref_dir=${bryce_data_dir}/prediction_files/gene_ref_files/chr${script_chr}
training_dir=${bryce_data_dir}/model_training/chr${script_chr}
varID_dir=${training_dir}/varIDs
output_dosage_dir=${bryce_data_dir}/prediction_files/subset_dosage/chr${script_chr}


if [ $( awk -v chunk=$script_chunk '$5==chunk' ${gencode_ref_dir}/gencode_build38_chr${script_chr}_genes.txt | wc -l ) -eq 0 ]; then
	echo There are no genes on chromosome $script_chr in chunk $script_chunk. Make sure this is correct via the reference genecode file.
	exit 0
fi


while read chr start_pos end_pos gene_name chunk
do
	# IF the model was fit with Elastic Net
	if [ -f ${training_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection.betas.EN.txt ]; then
		#Output a varFile for any variant which was included in a gene expression model in the chr*chunk pair. 
		tail -n +2 ${training_dir}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection.betas.EN.txt | cut -f 1-4 | sed 's/^chr//g' | sed 's/\t/:/g' >> ${varID_dir}/DGN_chr${chr}_chunk${chunk}_UKB_intersection_varIDs_b38
	fi
	
	
done < <( awk -v chunk=$script_chunk '$5==chunk' ${gencode_ref_dir}/gencode_build38_chr${script_chr}_genes.txt ) 	

if [ ! -f  ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b38 ]; then
	echo No genes from this chunk were successfully fit
	exit 0
fi

#Get unique IDs
sort ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b38 | uniq > ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b38_uniq

#Lift them over

~/bin/id_liftOver.sh ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b38_uniq ~/bin/chain_files/hg38ToHg19.over.chain ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_b37_uniq ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b37_uniq 

#cleanup
rm ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b38

#Create variant inclusion file for the dosage computation
paste ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b37_uniq <( awk -F ":" ' { print "NA",$1,$2,$3,$4 }' OFS=" " ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b37_uniq ) -d ' ' | sed '1 i\SNPID rsid chromosome position alleleA alleleB' > ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varFile_b37_uniq.txt

#correct formatting issue for UKB data
if [ "$script_chr" -lt 10 ]; then
	awk '!($1 == "SNPID"){ $1=0$1; $3=0$3} { print $0 } ' ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varFile_b37_uniq.txt > ${varID_dir}/temp_${script_chr}_${script_chunk}
	mv ${varID_dir}/temp_${script_chr}_${script_chunk} ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varFile_b37_uniq.txt
fi

#Create subset dosage file for UKB data for prediction 
qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${script_chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${script_chr}_v3_s487395.sample -compare-variants-by 'position,alleles' -incl-variants ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varFile_b37_uniq.txt -incl-samples ~/scr/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -og ${output_dosage_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_modelSNPs.dosage.gz

rm ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b37_uniq  ${varID_dir}/DGN_chr${script_chr}_chunk${script_chunk}_UKB_intersection_varIDs_b38_uniq 
