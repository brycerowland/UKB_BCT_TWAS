#!/bin/bash

module add qctool

GWAS_file=/proj/yunligrp/users/bryce/twas/UKB/data/BCX_GWAS_2020_known_variants/Vuckovic_SuppTables_clean.tsv
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/UKB_snpStats_reported_varlists

cat <( paste -d "\t" <( awk 'NR>1 {print $5,$6,$7,$8}' OFS=":"  $GWAS_file ) <(  awk 'NR>1 { print "NA",$5,$6,$7,$8 }' OFS="\t" $GWAS_file) ) \
	<( paste -d "\t" <( awk 'NR>1 {print $5,$6,$8,$7}' OFS=":"  $GWAS_file ) <(  awk 'NR>1 { print "NA",$5,$6,$8,$7 }' OFS="\t" $GWAS_file )  ) | \
	sort -k1V | uniq | awk '($3<10){ $1=0$1; $3=0$3 } {print $0}' | sed 's/\t/ /g' | sed '1 i\SNPID rsid chromosome position alleleA alleleB' >  ${out_dir}/BCX_reported_vars.txt

for chr in {1..22}; do
	awk -v chr=$chr 'NR ==1 || $3==chr' ${out_dir}/BCX_reported_vars.txt > ${out_dir}/BCX_reported_vars_chr${chr}.txt

	sbatch -J snp_stat_file_BCX_variants_chr${chr} -t 1:00:00 --mail-type END --mail-user bryce.rowland@unc.edu --out /proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/out_files/snp_stat_file_BCX_variants_chr${chr}.out --wrap "qctool -g /proj/yunligrp/ukbiobank/genetics/data/imputed_v3/data/ukb_imp_chr${chr}_v3.bgen -s /proj/yunligrp/UKBB_phen_29983/sample_files_new/ukb25953_imp_chr${chr}_v3_s487395.sample -compare-variants-by 'position,alleles' -incl-variants ${out_dir}/BCX_reported_vars_chr${chr}.txt -incl-samples ~/scr/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs -snp-stats -osnp ${out_dir}/UKB_SnpStats_chr${chr}_BCX_variants.txt.gz" 
	
done


