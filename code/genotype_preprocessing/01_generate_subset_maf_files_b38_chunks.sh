#!/bin/bash

data_dir="/pine/scr/b/r/bryce38/twas/UKB_temp/data"
b37_var_dir=${data_dir}/genotype_data/DGN_b37_var_lists
pheno_dir=${data_dir}/phenotype_data

module add qctool

#Note that in the UKB .bgen files, you need to write the chromosome as "0#" for chromosome #, when # < 10. But here this is already accounted for in the variant files, so only one
# loop is needed.
for chr in {22..1}; do
	for chunk in {42..1}; do
	for list_type in varFileFlipped varFile; do

	#Check to see if there is a corresponding chunk from the DGN eQTL files.	
	if [ -f ${b37_var_dir}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_${list_type}.txt ]
	then
		#If so, submit the job to compute snp stats on the UKB individuals with qctools. We only consider snps which are present in the eQTL dataset because of TWAS model training (we want the intersection of SNPs between the UKB data and DGN). 
		sbatch -J UKB_eur_snp_stats_chr${chr}_chunk${chunk}_${list_type} -t 6:00:00 --mem 1GB -o ../../data/imputed/out_files/UKB_eur_snp_stats_chr${chr}_chunk${chunk}_${list_type}.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample -compare-variants-by 'position,alleles' -incl-variants ${b37_var_dir}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_${list_type}.txt -incl-samples ${pheno_dir}/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -snp-stats -osnp ${data_dir}/genotype_data/maf_files/chr${chr}/chr${chr}_chunk${chunk}_snp_stats_intersectDGN_${list_type}"
	fi
  done
 done
done


#Rerun specific jobs.

: ' 
for chr in 2; do
        for chunk in 29; do
        for list_type in varFile; do

        #Check to see if there is a corresponding chunk from the DGN eQTL files.        
        if [ -f ${b37_var_dir}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_${list_type}.txt ]
        then
                #If so, submit the job to compute snp stats on the UKB individuals with qctools. We only consider snps which are present in the eQTL dataset because of TWAS model training (we want the intersection of SNPs between the UKB data and DGN). 
                sbatch -J UKB_eur_snp_stats_chr${chr}_chunk${chunk}_${list_type} -t 8:00:00 --mem 2GB -o ../../data/imputed/out_files/UKB_eur_snp_stats_chr${chr}_chunk${chunk}_${list_type}.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample -compare-variants-by 'position,alleles' -incl-variants ${b37_var_dir}/chr${chr}/DGN_b37_chr${chr}_chunk${chunk}_${list_type}.txt -incl-samples ${pheno_dir}/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -snp-stats -osnp ${data_dir}/genotype_data/maf_files/chr${chr}/chr${chr}_chunk${chunk}_snp_stats_intersectDGN_${list_type}"
        fi
  done
 done
done

'
