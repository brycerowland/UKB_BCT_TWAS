#!/bin/bash

UKB_pheno_dir="/proj/yunligrp/UKBB_phen_29983"
out_dir="/pine/scr/b/r/bryce38/twas/UKB_temp/data/phenotype_data"

#Select relevant columns
cut ${UKB_pheno_dir}/ukb32796.csv -f 1,2 -d ',' | sed 's/,/\t/g' > ${out_dir}/ukb32796_BCT_TWAS_cols.tsv
cut ${UKB_pheno_dir}/ukb28821.tab -f 1,23,749-750,758-797,904,917,943,930,969,982,956,1008,995,1034 > ${out_dir}/ukb28821_BCT_TWAS_cols.tsv
cut ${UKB_pheno_dir}/ukb26645.tab -f 1,59,65,68,71,74,77,80,86,89,92,95,98,104,107,110,113,116,119,122 > ${out_dir}/ukb26645_BCT_TWAS_cols.tsv
