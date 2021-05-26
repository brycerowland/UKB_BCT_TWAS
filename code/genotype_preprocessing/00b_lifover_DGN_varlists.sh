#!/bin/bash

out_path="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data/DGN_b38_var_lists"
b37_out_path="/pine/scr/b/r/bryce38/twas/UKB_temp/data/genotype_data/DGN_b37_var_lists"

for chr in {1..22}; do
	out_dir=${out_path}/chr${chr}	
	echo chr${chr}. . .
	for chunk in {1..42}; do
	if [ -f ${out_dir}/DGN_b38_chr${chr}_chunk${chunk}_varIDs ]	
	then		
		
		#Remove "chr"
		sed -i 's/^chr//g' ${out_dir}/DGN_b38_chr${chr}_chunk${chunk}_varIDs

		#Liftover IDs using Bryce"s script. 
		~/bin/id_liftOver.sh ${out_dir}/DGN_b38_chr${chr}_chunk${chunk}_varIDs  ~/bin/chain_files/hg38ToHg19.over.chain ${out_dir}/DGN_b38_chr${chr}_chunk${chunk} ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varIDs

		#Create variant inclusion file for the dosage computation
        	paste ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varIDs <( awk -F ":" ' { print "NA",$1,$2,$3,$4 }' OFS=" " ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varIDs ) -d ' ' | sed '1 i\SNPID rsid chromosome position alleleA alleleB' > ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varFile.txt

        	if [ "$chr" -lt 10 ]; then
                	awk '!($1 == "SNPID"){ $1=0$1; $3=0$3} { print $0 } ' ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varFile.txt > pizza
                	mv pizza ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varFile.txt
        	fi

		#Create flipped version. To check if variants are present in UKB, but alleles are flipped compared to DGN.
                awk -F ":" '{ print $1":"$2":"$4":"$3,"NA",$1,$2,$4,$3 }' ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varIDs | sed '1 i\SNPID rsid chromosome position alleleA alleleB' > ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varFileFlipped.txt

                if [ "$chr" -lt 10 ]; then
                        awk '!($1 == "SNPID"){ $1=0$1; $3=0$3} { print $0 } ' ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varFileFlipped.txt > pizza
                        mv pizza ${b37_out_path}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_varFileFlipped.txt
                fi

	fi

	done

done

