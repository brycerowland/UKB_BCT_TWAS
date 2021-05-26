#!/bin/bash

gene_ref_file=/proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/gene_ref_files/gencode_build38_allChr_TWAS_genes.txt
gexp_dir=/proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/predicted_gexp

chr=$1

rm -f  ${gexp_dir}/DGN_UKB_intersection_chr${chr}_trickDosage.gz

while read chr start_pos end_pos gene_name chunk en_fit; do 
	echo $gene_name
	midpoint=$( echo $(( ($start_pos + $end_pos)/2 )) )	

	#Assemble the predicted gene expression file into the 'dosage' format for imputed SNPs in boltLMM.
        paste <( echo "$gene_name $chr $midpoint A T" ) <( cut -f 2 -d ' ' ${gexp_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_${gene_name}_UKB_intersection.grex2 | tawk) -d ' ' | gzip -c >> ${gexp_dir}/DGN_UKB_intersection_chr${chr}_trickDosage.gz

done < <( tail -n +2 $gene_ref_file | awk -v chr=$chr '$1==chr' ) 
