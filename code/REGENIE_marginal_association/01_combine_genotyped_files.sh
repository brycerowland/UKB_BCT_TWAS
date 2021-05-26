#!/bin/bash

module add plink 
rm -f list_beds.txt

plink_dir=/proj/yunligrp/users/bryce/twas/UKB/data/minority_TWAS_replication/UKB_LD_pruned_genotypes

for chr in {2..22}; do echo "${plink_dir}/ukb_chr${chr}_LD_pruned_MAF0.01.bed ${plink_dir}/ukb_chr${chr}_LD_pruned_MAF0.01.bim ${plink_dir}/ukb_chr${chr}_LD_pruned_MAF0.01.fam" >> list_beds.txt; done

plink \
  --bfile ${plink_dir}/ukb_chr1_LD_pruned_MAF0.01 \
  --merge-list list_beds.txt \
  --make-bed --out ${plink_dir}/ukb_allChrs_LD_pruned_MAF0.01

