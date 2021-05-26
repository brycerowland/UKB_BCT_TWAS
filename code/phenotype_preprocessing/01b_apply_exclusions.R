#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

write_id_file <- function(data, name){
	data %>%
        as.data.frame() %>%
        write_tsv(paste0("/nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data/", name), col_names = F)
}

nih_dir <- "/nas/longleaf/home/bryce38/bryce_group/twas/UKB/data/phenotype_data/nih_files"
pheno_dir_scr <- "/nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data"

#Files generated from bash one liners, subsets of ub32796.csv
drugs <- fread(paste(pheno_dir_scr, "ukb32796_meds.csv", sep = "/"), header = T, data.table = F)
c20001 <- fread("/nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data/ukb32796_cancer.csv", header = T, data.table = F)
c41202 <- fread("/nas/longleaf/home/bryce38/scr/twas/UKB_temp/data/phenotype_data/ukb32796_icd10_main.csv", header = T, data.table = F)
c41204 <- fread(paste(pheno_dir_scr, "ukb32796_icd10_secondary.csv", sep = "/"), header = T, data.table = F)
c41203 <- fread(paste(pheno_dir_scr, "ukb32796_icd9_main.csv", sep = "/"), header = T, data.table = F)
c41205 <- fread(paste(pheno_dir_scr, "ukb32796_icd9_secondary.csv", sep = "/"), header = T, data.table = F)
c41200 <- fread(paste(pheno_dir_scr, "ukb32796_surgical_main.csv", sep = "/"), header = T, data.table = F)
c41210 <- fread(paste(pheno_dir_scr, "ukb32796_surgical_secondary.csv", sep = "/"), header = T, data.table = F)

#Files from NIH
code19 <- fread(paste(nih_dir, "coding19.tsv", sep = "/"),header=T,data.table=F)
code240 <- fread(paste(nih_dir, "coding240.tsv", sep = "/"),header=T,data.table=F)
code87 <- fread(paste(nih_dir, "coding87.tsv", sep = "/"),header=T,data.table=F)

#Get ID's for exclusion based on drugs
drugs.excl <- apply(drugs[,-1],1,function(x)any(x%in%c(1140869523:1140870285,1141157448,1140870568,1140870600,1140870604,1140870618,1140870626,1141178816,1141178858),na.rm=T))
drugs.excl <- drugs[drugs.excl,1]

#Get ID's for exclusion based on self reported cancer
ccode <- c(1085,1074,1070,1063,1058,1056,1055,1053,1052,1051,1050,1048,1047)
c20001.excl <- apply(c20001[,-1],1,function(x)any(x%in%ccode,na.rm=T))
c20001.excl <- c20001[c20001.excl,1]

##ICD10
#Convert codes (provided by NIH)
icd10.main <- c("B20","B21","B22","B23","B24","C40","C41","C81","C82","C83","C84","C85","C86","C87","C88","C89","C90","C91","C92","C93","C94","C95","C96","D45","D46","D47","D55","D56","D57","D58","D59","D60","D61","D63",
              "D640","D641","D642","D643","D644","D65","D66","D67","D68","D69","D70","D71","D72","D73","D74","D75","D76","D77","D80","D81","D82","D83","D84","D85","D86","D87","D88","D89","K70","K71","K74",
              "R70","R71","R72","R73","R74","R75","R76","R77","R78","R79")
icd10 <- NULL

for (i in icd10.main) icd10 <- c(icd10,code19[grep(i,code19[,1]),1])
icd10 <- icd10[-grep("Block",icd10)]

#Get ID's for exclusion based on ICD10 - main
c41202.excl <- apply(c41202[,-1],1,function(x)any(x%in%icd10,na.rm=T))
c41202.excl <- c41202[c41202.excl,1]

#Get ID's for exclusion based on ICD10 - secondary
c41204.excl <- apply(c41204[,-1],1,function(x)any(x%in%icd10,na.rm=T))
c41204.excl <- c41204[c41204.excl,1]

##ICD9
#Convert codes (provided by NIH)
icd9.main <- c(170,200:208,2384,2385,2386,2387,286,287,288,289,282,283,284,571,790)
icd9 <- code87[c(1205:1215,1403:1460,1675:1678,2194:2259,2271:2337,3958:3972,7540:7550),1]

#Get ID's for exclusion based on ICD9 - main 
c41203.excl <- apply(c41203[,-1],1,function(x)any(x%in%icd9,na.rm=T))
c41203.excl <- c41203[c41203.excl,1]

#Get ID's for exclusion based on ICD9 - secondary 
c41205.excl <- apply(c41205[,-1],1,function(x)any(x%in%icd9,na.rm=T))
c41205.excl <- c41205[c41205.excl,1]

##Operative Procedures
#convert codes (provdied by NIH)
op <- code240[c(3605:3623,8434:8442),1]

#Get ID's for exclusion based on Operative Procedures - main
c41200.excl <- apply(c41200[,-1],1,function(x)any(x%in%op,na.rm=T))
c41200.excl <- c41200[c41200.excl,1]

#Get ID's for exclusion based on Operative Procedures - secondary
c41210.excl <- apply(c41210[,-1],1,function(x)any(x%in%op,na.rm=T))
c41210.excl <- c41210[c41210.excl,1]

#Write ID files
write_id_file(drugs.excl, "drugs_exclusions.id")
write_id_file(c20001.excl, "cancer_exclusions.id")
write_id_file(c41202.excl, "icd10_main_exclusions.id")
write_id_file(c41204.excl, "icd10_secondary_exclusions.id")
write_id_file(c41203.excl, "icd9_main_exclusions.id")
write_id_file(c41205.excl, "icd9_secondary_exclusions.id")
write_id_file(c41200.excl, "operative_procedures_main_exclusions.id")
write_id_file(c41210.excl, "operative_procedures_secondary_exclusions.id")

