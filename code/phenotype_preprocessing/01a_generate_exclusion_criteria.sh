#!/bin/bash

module add r/3.6.0

#Purpose of this file is to:
#	1. Demonstrate where each filtering criteria was found in the columns of the phenotype file
#	2. Generate a list of Europeans who past exclusion criteria in UKB
# This file calls a seperate R script which reads in some of the subset files to perform efficient searching with `apply`.

pheno_file="/proj/yunligrp/UKBB_phen_29983/ukb32796.csv"
pheno_data_dir="/nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data"

#Pregnancy - column 359
: '
grep -w "3140" ukb32796_cols 
   359	3140-0.0
   360	3140-1.0
   361	3140-2.0
  3140	20229-0.39
'
awk -F "," '$359 == "\"1\"" { print $1 }' ${pheno_file} | sed 's/"//g' > ${pheno_data_dir}/preg_exclusions.id


#Cancer treatments/drugs
: ' 
#First, search for medication columns at baseline.
grep -w "20003-0" ukb32796_cols
   902	20003-0.0
   903	20003-0.1
   904	20003-0.2
   905	20003-0.3
   906	20003-0.4
   907	20003-0.5
   943	20003-0.41
...
   944	20003-0.42
   945	20003-0.43
   946	20003-0.44
   947	20003-0.45
   948	20003-0.46
   949	20003-0.47
' 
#Create subset file for later R processing
cut -f 1,902-949 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_meds.csv

: '
grep -w "20001-0" ukb32796_cols 
   782	20001-0.0
   783	20001-0.1
   784	20001-0.2
   785	20001-0.3
   786	20001-0.4
   787	20001-0.5
'

#Self-reported Cancer Status
cut -f 1,782-787 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_cancer.csv

#ICD10 main and secondary diag
: '
grep -w "41202" ukb32796_cols 
  5685	41202-0.0
  5686	41202-0.1
  5687	41202-0.2
  5688	41202-0.3
  5689	41202-0.4
  5745	41202-0.60
...
  5746	41202-0.61
  5747	41202-0.62
  5748	41202-0.63
  5749	41202-0.64
  5750	41202-0.65
'

: ' 
grep -w "41204" ukb32796_cols 
  5779	41204-0.0
  5780	41204-0.1
  5781	41204-0.2
  5782	41204-0.3
  5783	41204-0.4
  5784	41204-0.5
...
  5955	41204-0.176
  5956	41204-0.177
  5957	41204-0.178
  5958	41204-0.179
  5959	41204-0.180
  5960	41204-0.181
  5961	41204-0.182
  5962	41204-0.183
'

cut -f 1,5685-5750 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_icd10_main.csv
cut -f 1,5779-5962 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_icd10_secondary.csv

#ICD-9 main and secondary diagnoses

: '
grep -w "41203" ukb32796_cols 
  5751	41203-0.0
  5752	41203-0.1
  5753	41203-0.2
  5754	41203-0.3
  5755	41203-0.4
  5756	41203-0.5
...  
  5772	41203-0.21
  5773	41203-0.22
  5774	41203-0.23
  5775	41203-0.24
  5776	41203-0.25
  5777	41203-0.26
  5778	41203-0.27
'

: ' 
grep -w "41205" ukb32796_cols 
  5963	41205-0.0
  5964	41205-0.1
  5965	41205-0.2
  5966	41205-0.3
  5967	41205-0.4
  5968	41205-0.5
...
  5987	41205-0.24
  5988	41205-0.25
  5989	41205-0.26
  5990	41205-0.27
  5991	41205-0.28
  5992	41205-0.29
'

cut -f 1,5751-5778 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_icd9_main.csv
cut -f 1,5963-5992 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_icd9_secondary.csv

#Surgical Procedures - main and secondary

: ' 
grep -w "41200" ukb32796_cols 
  5614	41200-0.0
  5615	41200-0.1
  5616	41200-0.2
  5617	41200-0.3
  5618	41200-0.4
  5619	41200-0.5
  5620	41200-0.6 
...
  5661	41200-0.47
  5662	41200-0.48
'

: ' 
grep -w "41210" ukb32796_cols 
  6048	41210-0.0
  6049	41210-0.1
  6050	41210-0.2
  6051	41210-0.3
  6052	41210-0.4
  6053	41210-0.5
...
  6129	41210-0.81
  6130	41210-0.82
  6131	41210-0.83
  6132	41210-0.84
  6133	41210-0.85
'

cut -f 1,5614-5662 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_surgical_main.csv
cut -f 1,6048-6133 -d "," ${pheno_file} > ${pheno_data_dir}/ukb32796_surgical_secondary.csv

#European IDs
#awk -F "," '$737 == "\"1003\"" || $737 == "\"1001\"" || $737 == "\"1002\"" || $737 == "\"1\"" {print $1"\t"$737}'  /proj/yunligrp/UKBB_phen_29983/ukb28821.csv  | sed 's/"//g' > ${pheno_data_dir}/ukb_euro_ids_w_ethnicity_code

#Rather than doing European ancestry based off of self report - we use the k-means+Self Report Exclusions clustering approach as determined by Misa Graff. We also filter out any mixed ancestry individuals from self-report data. 
awk -F "\t" '$67 == 1 && $15!~/^2.*/ { print $1 }' /proj/yunligrp/UKBB_phen_29983/example_scripts_GRM_Misa/covariate_file/UKBb_allcovar_geneticsample_12Feb2020.txt > ${pheno_data_dir}/ukb_euro_ids_kMeans_selfReport

#poor heterozygosity/missingness (Recommended genomic analysis exclusions, and genetic relatedness exclusions)
tail -n +2 /proj/yunligrp/UKBB_phen_29983/ukb28821.tab | awk -F "\t" '$814 == 1 || $814 ==2 || $798 == 1 { print $1}' > ${pheno_data_dir}/genome_exclusions.id

#Withdrawn consent
cat /proj/yunligrp/UKBB_phen_29983/samples_to_exclude/* | sort | uniq > ${pheno_data_dir}/withdrawn_concent_exclusions.id

#Run R script to process new .csv
echo Running R script to process exclusions
./01b_apply_exclusions.R

#Get unique exclusions
cat ${pheno_data_dir}/*_exclusions.id | sort | uniq > ${pheno_data_dir}/unique_exclusions.id

#Apply exclusions
awk 'NR==FNR {a[$1]; next} !($1 in a)' ${pheno_data_dir}/unique_exclusions.id ${pheno_data_dir}/ukb_euro_ids_kMeans_selfReport  > ${pheno_data_dir}/ukb_euro_ids_kMeans_selfReport_ExclusionsApplied

