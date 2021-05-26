#!/bin/bash

module add qctool
#Note that in the UKB .bgen files, you need to write the chromosome as "0#" for chromosome #, when # < 10.
for chr in {1..9}; do
	for i in {1..42}; do
	#sbatch -J UKB_eur_dosage_${chr}_chunk${i} -t 3-00:00:00 --mem 8GB -o ../../data/imputed/out_files/UKB_eur_dosage_chr${chr}_chunk${i}.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample -incl-range 0${chr}:$(( ($i-1)*6000000 + 1 ))-$(( $i*6000000 + 4000000 )) -og /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/genotype_data/eur_dosage/chr${chr}/UKB_imp_chr${chr}_chunk${i}.dosage.gz -incl-samples /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs" 
	sbatch -J UKB_eur_maf_chr${chr}_chunk${i}_PiHat0.2 -t 1-00:00:00 --mem 8GB -o ../../data/imputed/out_files/UKB_eur_maf_chr${chr}_chunk${i}_PiHat0.2.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample  -incl-range 0${chr}:$(( ($i-1)*6000000 + 1 ))-$(( $i*6000000 + 4000000 )) -incl-samples /pine/scr/b/r/bryce38/twas/UKB_temp/data/GRM/UKB_BCT_GRM_PiHat0.2_NEW_IDs -snp-stats -osnp /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/genotype_data/PiHat0.2/maf_files/chr${chr}/chr${chr}_chunk${i}_snp_stats_PiHat0.2"
 done
done

for chr in {10..22}; do
        for i in {1..42}; do
	#sbatch -J UKB_eur_dosage_${chr}_chunk${i} -t 3-00:00:00 --mem 8GB -o ../../data/imputed/out_files/UKB_eur_dosage_chr${chr}_chunk${i}.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample -incl-range ${chr}:$(( ($i-1)*6000000 + 1 ))-$(( $i*6000000 + 4000000 )) -og /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/genotype_data/eur_dosage/chr${chr}/UKB_imp_chr${chr}_chunk${i}.dosage.gz -incl-samples /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs"
	sbatch -J UKB_eur_maf_chr${chr}_chunk${i}_PiHat0.2 -t 1-00:00:00 --mem 8GB -o ../../data/imputed/out_files/UKB_eur_maf_chr${chr}_chunk${i}_PiHat0.2.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample  -incl-range ${chr}:$(( ($i-1)*6000000 + 1 ))-$(( $i*6000000 + 4000000 )) -incl-samples /pine/scr/b/r/bryce38/twas/UKB_temp/data/GRM/UKB_BCT_GRM_PiHat0.2_NEW_IDs -snp-stats -osnp /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/genotype_data/PiHat0.2/maf_files/chr${chr}/chr${chr}_chunk${i}_snp_stats_PiHat0.2"
 done
done

: '
for chr in 14; do
        for i in 11; do
        sbatch -J UKB_eur_dosage_${chr}_chunk${i}_T -t 1-00:00:00 --mem 8GB -o ../../data/TNFAIP2/new_dosage/out_files/UKB_eur_dosage_chr${chr}_chunk${i}_T.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample -incl-range ${chr}:102000000-105000000 -og /proj/yunligrp/users/bryce/twas/UKB/data/TNFAIP2/new_dosage/UKB_imp_chr${chr}_chunk${i}.dosage.gz -incl-samples /pine/scr/b/r/bryce38/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs"
	sbatch -J UKB_eur_maf_chr${chr}_chunk${i}_T -t 8:00:00 --mem 8GB -o ../../data/TNFAIP2/new_dosage/out_files/UKB_eur_maf_chr${chr}_chunk${i}_T.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample  -incl-range ${chr}:102000000-105000000 -incl-samples /pine/scr/b/r/bryce38/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -snp-stats -osnp /proj/yunligrp/users/bryce/twas/UKB/data/TNFAIP2/new_dosage/chr${chr}_chunk${i}_snp_stats"
 done
done
'
: '
for chr in 16; do
        for i in 9; do
        sbatch -J UKB_eur_maf_chr${chr}_chunk${i} -t 5-00:00:00 --mem 8GB -o ../../data/imputed/out_files/UKB_eur_maf_chr${chr}_chunk${i}_test.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample  -incl-range ${chr}:$(( ($i-1)*10000000 + 1 ))-$(( $i*10000000 )) -incl-samples /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -snp-stats -osnp /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/genotype_data/maf_files/chr${chr}/chr${chr}_chunk${i}_snp_stats"
 done
done
'
: '
while read -r chr i; do
	sbatch -J UKB_eur_maf_chr${chr}_chunk${i} -t 2-00:00:00 --mem 8GB -o ../../data/imputed/out_files/UKB_eur_maf_chr${chr}_chunk${i}_test.out --mail-type END --mail-user bryce.rowland@unc.edu --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample  -incl-range ${chr}:$(( ($i-1)*10000000 + 1 ))-$(( $i*10000000 )) -incl-samples /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -snp-stats -osnp /nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/genotype_data/maf_files/chr${chr}/chr${chr}_chunk${i}_snp_stats"
done <  <( tail -n +2 maf_jobs_to_be_rerun | awk '$1 >= 10' )
'
