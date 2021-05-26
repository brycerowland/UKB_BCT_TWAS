#!/bin/bash

#The purpose of this file is to get variants which are well imputed (MAF > 0.05 and R2 > .8) in UKB, and then get this list to subset the imputed DGN data. 
# We output three lists of variants, non-flipped variant ID's, flipped variant ID's, and the joint list. We output seperate flipped and non-flipped lists for going
# back to UKB after model training. 

data_dir="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data"
snp_dir=${data_dir}/maf_files
liftover_input_dir=${data_dir}/UKB_b37_liftover_input
liftover_output_dir=${data_dir}/UKB_b38_liftover_output
DGN_dir="/proj/yunligrp/users/jwen/copy_scr_munan/DGN/cohort_QC/DGN/chunking.raw"
output_dir=${data_dir}/DGN_UKB_filtered_imputed_data/DGN_filtered_vcf

module add samtools

chr=$1

for chunk in {42..1}; do
	#First, check both for flipped and non-flipped allele files. If they exist, print out variant IDs in DGN's orientation.

	# If the non-flipped snp stats file exists
	if [ -f ${snp_dir}/chr${chr}/chr${chr}_chunk${chunk}_snp_stats_intersectDGN_varFile ]; then
		#Field 14 is minor allele frequency, and is always the minimum between field 12 and 13 (allele A and allele B frequency). 
		grep -v "^#" ${snp_dir}/chr${chr}/chr${chr}_chunk${chunk}_snp_stats_intersectDGN_varFile | tail -n +2 |awk -v chr=$chr '$14 > 0.05 && $18 > .8 { print chr":"$4":"$5":"$6 }' > ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_nonFlipped
	fi


	# If the flipped snp stats file exists, we output an allele ID file with alleles in the orientation of DGN. That is to say, we flip them back. 
	if [ -f ${snp_dir}/chr${chr}/chr${chr}_chunk${chunk}_snp_stats_intersectDGN_varFileFlipped ]; then
		#Field 14 is minor allele frequency, and is always the minimum between field 12 and 13 (allele A and allele B frequency). 
                        grep -v "^#" ${snp_dir}/chr${chr}/chr${chr}_chunk${chunk}_snp_stats_intersectDGN_varFileFlipped | tail -n +2 | awk -v chr=$chr '$14 > 0.05 && $18 > .8 { print chr":"$4":"$6":"$5 }' > ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_Flipped	
	fi

	#Join them
	if [ -f ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_Flipped ] || [ -f ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_Flipped ]; then
		#This find command checks to see if there are files which match this pattern in the given directory, and then uses the exec command to cat them together. 
		# This handles the situation where only one snp_stat file may be generated.
		find ${liftover_input_dir}/chr${chr} -name "DGN_chr${chr}_chunk${chunk}_*Flipped" -exec cat {} \;>  ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_wMultiAlleles

		#Remove multiallelic sites in UKB
		awk -F ":" 'NR==FNR {a[$2]++; next} a[$2] < 2' ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_wMultiAlleles ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_wMultiAlleles > ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps

		#Liftover using the ID lift script
                	~/bin/id_liftOver.sh ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps ~/bin/chain_files/hg19ToHg38.over.chain ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38snps ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs

		#Add "chr" back in
		sed -i 's/^/chr/g' ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs

		#Go back and remove multiallelic variants from the Flipped and nonFlipped files if there are any. This is for later use of these files.
		if [ $( cat ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_wMultiAlleles | wc -l) -ne $( cat ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps  | wc -l ) ]; then
			echo $chr $chunk

			#Fix flipped
			awk 'NR==FNR {a[$1];next} $1 in a' ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_Flipped > ${liftover_input_dir}/chr${chr}/temp
			mv ${liftover_input_dir}/chr${chr}/temp ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_Flipped

			#Fix non-flipped	
			awk 'NR==FNR {a[$1];next} $1 in a' ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_nonFlipped > ${liftover_input_dir}/chr${chr}/temp
                                mv ${liftover_input_dir}/chr${chr}/temp ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_nonFlipped

		fi

	       #Now, we lift over the individuals var lists to b38. 
                       #Lift it over
                       ~/bin/id_liftOver.sh ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_nonFlipped ~/bin/chain_files/hg19ToHg38.over.chain ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38snps_nonFlipped ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_nonFlipped
 
                       #Add "chr" back in
                       sed -i 's/^/chr/g' ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_nonFlipped

                       #lift them over
                        ~/bin/id_liftOver.sh ${liftover_input_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b37snps_Flipped ~/bin/chain_files/hg19ToHg38.over.chain ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38snps_Flipped ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_Flipped
 
                       #Add "chr" back in
                        sed -i 's/^/chr/g' ${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs_Flipped


		#Subset the vcf file on the b38 ID's
                bcftools filter -i "ID=@${liftover_output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_b38_IDs" -Oz -o ${output_dir}/chr${chr}/DGN_chr${chr}_chunk${chunk}_UKB_intersection.vcf.gz ${DGN_dir}/chunking.chr${chr}/DGN.chr${chr}.chunk${chunk}.vcf.gz

	fi
done
