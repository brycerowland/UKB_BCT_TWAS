#!/bin/bash

data_dir=/proj/yunligrp/users/bryce/twas/UKB/data
gexp_dir=${data_dir}/prediction_files/predicted_gexp
out_dir=${data_dir}/conditional_analysis/REGENIE/sig_gene_predicted_gexp

module add plink/2.00a-20190527

pheno_class=$1
chr=$2
 
sig_gene_list=${data_dir}/analysis_results/regenie/association_results/${pheno_class}_TWAS_bonferroniSigGenes
awk 'NR==FNR {a[$4]; next} $1 in a' <( awk -v chr=$chr '$1==chr' $sig_gene_list ) <( zcat ${gexp_dir}/DGN_UKB_intersection_chr${chr}_trickDosage.gz ) > ${out_dir}/${pheno_class}_TWAS_bonferroniSigGenes_predGexp_chr${chr}

plink2 --import-dosage ${out_dir}/${pheno_class}_TWAS_bonferroniSigGenes_predGexp_chr${chr} noheader chr-col-num=2 pos-col-num=3 skip1=2 --fam ${gexp_dir}/DGN_UKB_intersection_trickDosage_IDs.fam --hard-call-threshold 0 --out ${out_dir}/${pheno_class}_TWAS_bonferroniSigGenes_predGexp_chr${chr}

