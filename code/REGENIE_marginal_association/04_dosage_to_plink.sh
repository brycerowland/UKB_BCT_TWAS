#!/bin/bash

#This file converts the TrickDosage files (which were generated for use in bolt-LMM but is basically plink1 dosage) to plink2 files which store the dosage information. These files are compatible with REGENIE for the analysis. 

#Likely, if we go this route in the future, I will need a script which goes from the raw predicted gene expression values to the plink2 format, but  this is using what we currently have on hand. 

module add plink/2.00a-20190527 # has to be this version of plink

gexp_dir=/proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/predicted_gexp
: ' 
for chr in {1..22}; do
	zcat ${gexp_dir}/DGN_UKB_intersection_chr${chr}_trickDosage.gz
done > ${gexp_dir}/DGN_UKB_intersection_trickDosage
'
plink2 --import-dosage ${gexp_dir}/DGN_UKB_intersection_trickDosage noheader chr-col-num=2 pos-col-num=3 skip1=2 --fam ${gexp_dir}/DGN_UKB_intersection_trickDosage_IDs.fam --hard-call-threshold 0 --out ${gexp_dir}/DGN_UKB_intersection_allChr_trickDosage

#rm ${gexp_dir}/DGN_UKB_intersection_trickDosage
