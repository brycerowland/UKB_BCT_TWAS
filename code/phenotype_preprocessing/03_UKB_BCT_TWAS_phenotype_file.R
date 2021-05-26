#!/usr/bin/env Rscript

library(tidyverse)

theme_set(theme_bw())

pheno_dir <- "/pine/scr/b/r/bryce38/twas/UKB_temp/data/phenotype_data/"
out_dir <- "/proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/"


#Get european ids
eur_ids <- read_delim(str_c(pheno_dir, "ukb_euro_ids_kMeans_selfReport_ExclusionsApplied"),
                      col_names = F,
                      delim = ' ') %>% 
  select(id = X1)

sex <- read_tsv(str_c(pheno_dir, "ukb32796_BCT_TWAS_cols.tsv")) %>% 
  filter(eid %in% eur_ids$id)

cols <- read_tsv(str_c(pheno_dir, "ukb28821_BCT_TWAS_cols.tsv")) %>% 
  filter(f.eid %in% eur_ids$id)

wbc_extra <- read_tsv(str_c(pheno_dir, "ukb26645_BCT_TWAS_cols.tsv")) %>% 
  filter(f.eid %in% eur_ids$id)

#Vector for naming genotype variables
pc_names <- paste0("genotype_pc", 1:40)

pheno <- left_join(cols, sex, by = c("f.eid" = "eid")) %>% 
  left_join(wbc_extra, by = "f.eid") %>% 
  filter(complete.cases(.)) %>%
  #(n = 399849 after missing data)
  rename(sex = `31-0.0`, 
         age_at_recruitment = `f.21022.0.0`,
         white_blood_cell_count = `f.30000.0.0`,
         red_blood_cell_count = `f.30010.0.0`,
         hematocrit_percentage = `f.30030.0.0`,
         hemoglobin_concentration = `f.30020.0.0`,
         mean_cell_hemoglobin = `f.30050.0.0`,
         mean_cell_hemoglobin_concentration = `f.30060.0.0`,
         mean_cell_volume = `f.30040.0.0`,
         platelet_count = `f.30080.0.0`,
         red_blood_cell_width = `f.30070.0.0`,
         mean_platelet_volume = `f.30100.0.0`,
         genotype_measurement_batch = f.22000.0.0,
         center = f.54.0.0,
	 lymphocyte_count = f.30120.0.0,
	 monocyte_count = f.30130.0.0,
         basophil_count = f.30160.0.0,
         eosinophil_count = f.30150.0.0,
         neutrophil_count = f.30140.0.0,
	 platelet_distribution_width = f.30110.0.0,
	 plateletcrit = f.30090.0.0,
	 mean_sphered_cell_volume = f.30270.0.0,
	 reticulocyte_count = f.30250.0.0,
	 reticulocyte_percentage = f.30240.0.0,
	 immature_reticulocyte_fraction = f.30280.0.0,
	 hls_reticulocyte_count = f.30300.0.0,
	 hls_reticulocyte_percentage = f.30290.0.0,
	 mean_reticulocyte_volume = f.30260.0.0,
	 monocyte_percentage = f.30190.0.0,
	 neutrophil_percentage = f.30200.0.0,
	 eosinophil_percentage = f.30210.0.0,
	 basophil_percentage = f.30220.0.0,
	 lymphocyte_percentage = f.30180.0.0) %>% 
  filter(!(white_blood_cell_count  > 200 |
	platelet_count > 1000 |
	hemoglobin_concentration > 20 |
	hematocrit_percentage > 60 )) %>%
  mutate(eid = as.factor(f.eid))

colnames(pheno)[5:44] <- pc_names

pheno_ids <- pheno$f.eid %>% 
  as.data.frame()
nrow(pheno_ids)

#write files
write_tsv(pheno_ids, str_c(out_dir, "UKB_BloodCellTraits_EUR_TWAS_Analysis_IDs"),
          col_names = F)
write_tsv(pheno, str_c(out_dir, "UKB_BloodCellTraits_EUR_TWAS_w_varNames.tsv"))
write_tsv(pheno %>% filter(complete.cases(.)), str_c(out_dir, "UKB_BloodCellTraits_EUR_TWAS_w_varNames_CompleteCase.tsv"))
