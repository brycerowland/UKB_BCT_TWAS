library(tidyverse)


phenos <- read_tsv('../../../data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_w_varNames.tsv')


phenos %>% 
  select(-c(f.eid:genotype_pc40), -eid, -sex) %>% 
  summarise(across(.fns = list(mean = mean, 
                               sd = sd),
                   .names = "{.fn}_{.col}")) %>% 
  pivot_longer(cols = everything(), 
               names_to = c("measure", "phenotype"), 
               names_pattern = "([^_]*)_(.*)", 
               values_to = "value") %>% 
  pivot_wider(id_cols = phenotype, 
              names_from = measure, 
              values_from = value)%>% 
  transmute(Phenotype = str_replace_all(phenotype, 
                                    "_", " ") %>% 
           str_to_title(), 
         `Mean (SD)` = str_c(round(mean,3), " (", round(sd,3), ")")) %>% 
  write_tsv("../../../data/phenotype_data/SuppTable_Phenotype_distribution.tsv")
