#!/usr/bin/env Rscript

library(tidyverse)

#Define helper functions

inverse_normal_transformation <- function(r){
  qnorm((rank(r,na.last="keep")-0.5)/sum(!is.na(r)))
}

pheno_transform <- function(phenotype, .data, .log_list){
  log <- phenotype %in% .log_list
  if(log){
    str <- str_c("log10(", phenotype, " + 1) ~ " ,control_vars)
  } else {
    str <- paste(phenotype, control_vars, sep = "~")
  }
  f <- formula(str)
  r <- inverse_normal_transformation(resid(lm(f, data = .data, na.action = "na.exclude")))
  return(r)
}


#Read in data

pheno <- read_tsv("/proj/yunligrp/users/bryce/prs/trans_ethnic_PRS/data/BCT_phenotypes/UKB_BloodCellTraits_SAS_TWAS_w_varNames.tsv")


test_data <- pheno %>% 
  filter(!complete.cases(.))
#Create formula string for adjustment
control_vars <- (paste(c(paste0("genotype_pc", 1:10), 
                         c("sex", "age_at_recruitment", 
                           "I(age_at_recruitment^2)", "center",
                           "genotype_measurement_batch")), 
                       collapse = "+"))

# run analysis

pheno_list <- colnames(pheno)[c(45:54,56:74)]
log_list <- c("white_blood_cell_count", "lymphocyte_count", 
              "monocyte_count", "neutrophil_count", "eosinophil_count",
              "basophil_count", "lymphocyte_percentage", "monocyte_percentage",
              "neutrophil_percentage", "eosinophil_percentage", "basophil_percentage")
names(pheno_list) <- pheno_list

pheno_list

adjusted_BCT_phenotypes <- bind_cols(pheno %>% select(-any_of(pheno_list)), 
                                     map_dfc(pheno_list, ~pheno_transform(.x, pheno, log_list))) %>% 
  select(ID = f.eid, pheno_list)
nrow(adjusted_BCT_phenotypes)
nrow(pheno)
print("writing phenotype files")
write_tsv(adjusted_BCT_phenotypes, "/proj/yunligrp/users/bryce/prs/trans_ethnic_PRS/data/BCT_phenotypes/UKB_BloodCellTraits_SAS_TWAS_w_varNames_adjusted.tsv")
