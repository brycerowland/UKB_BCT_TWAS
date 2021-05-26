library(tidyverse)
library(UpSetR)
library(fastDummies)

ca_results <- read_tsv("../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_CA_results.tsv") %>% 
  filter(log10_conditional_p > -log10(0.05/11759)) %>% 
  select(gene_name, 
         phenotype)
fromList(ca_results)



p_dat <- ca_results %>% 
  mutate(phenotype = if_else(phenotype == "red_blood_cell_width", 
                             "Red_blood_cell_Distribution_Width", 
                             phenotype)) %>% 
  dummy_cols(select_columns = "phenotype") %>% 
  select(-phenotype) %>% 
  rename_with(.cols = starts_with("phenotype_"), 
              .fn = ~str_to_title(str_replace_all(str_remove(., "phenotype_"), "_", " "))) %>% 
  group_by(gene_name) %>% 
  summarise(across(.fns = sum)) %>% 
  as.data.frame()

upset(p_dat, 
      nsets = 10, 
      order.by = "freq")

