library(tidyverse)
library(UpSetR)
coloc_results <- read_tsv("../../../data/soranzo_group_results/COLOCS_Plus_Soranzo_annotation_Plus_Open_Targets.tsv") %>% 
  rename(twas_trait = phenotype)

all_pairs <- coloc_results %>% 
  distinct(phenotype_DEF, VAR)


coloc_assignments <- read_tsv("../gwas_variants_coloc_gene_assignments.tsv") %>% 
  rename(VAR = condsig_var, 
         twas_gene = gene_name, 
         twas_trait = UKB_TWAS_trait)

coloc_results %>% 
  distinct(VAR, twas_trait, COLOCS_gene)


coloc_csq <- coloc_results %>% 
  select(COLOCS_gene, VAR, twas_trait, Classif_CSQ, 
         soranzo_CSQ_gene, soranzo_GENE_EXP_gene) %>% 
  filter(soranzo_CSQ_gene != "NO_CANDIDATE") %>% 
  group_by(VAR, twas_trait)%>% 
  summarise(CSQ_at_least_one_match = sum(COLOCS_gene == soranzo_CSQ_gene)) %>% 
  ungroup() %>% 
  count(CSQ_at_least_one_match > 0) %>% 
  mutate(pct = n/sum(n)) %>% 
  slice(2) 

csq_pairs <-coloc_results %>% 
  select(COLOCS_gene, VAR, twas_trait, Classif_CSQ, 
         soranzo_CSQ_gene, soranzo_GENE_EXP_gene) %>% 
  mutate(soranzo_CSQ_gene= if_else(soranzo_CSQ_gene == "NO_CANDIDATE", 
                                   NA_character_, soranzo_CSQ_gene) ) %>% 
  
  group_by(VAR, twas_trait)%>% 
  summarise(CSQ_at_least_one_match = sum(COLOCS_gene == soranzo_CSQ_gene)) %>% 
  mutate(csq = if_else(CSQ_at_least_one_match > 0, 1, 0)) %>% 
  select(-CSQ_at_least_one_match)


coloc_gene_exp <- coloc_results %>% 
  select(COLOCS_gene, VAR, twas_trait, Classif_CSQ, 
         soranzo_CSQ_gene, soranzo_GENE_EXP_gene) %>% 
  filter(soranzo_GENE_EXP_gene != "NO_CANDIDATE") %>% 
  group_by(VAR, twas_trait) %>% 
  summarise(GENE_EXP_at_least_one_match = sum(COLOCS_gene == soranzo_GENE_EXP_gene)) %>% 
  ungroup() %>% 
  count(GENE_EXP_at_least_one_match > 0) %>% 
  mutate(pct = n/sum(n)) %>% 
  slice(2) 

gene_exp_pairs <- coloc_results %>%
  select(COLOCS_gene, VAR, twas_trait, Classif_CSQ, 
         soranzo_CSQ_gene, soranzo_GENE_EXP_gene) %>% 
  mutate(soranzo_GENE_EXP_gene= if_else(soranzo_GENE_EXP_gene == "NO_CANDIDATE", 
                                        NA_character_, soranzo_GENE_EXP_gene) ) %>%
  group_by(VAR, twas_trait) %>% 
  summarise(GENE_EXP_at_least_one_match = sum(COLOCS_gene == soranzo_GENE_EXP_gene)) %>% 
  ungroup() %>% 
  mutate(gexp = if_else(GENE_EXP_at_least_one_match > 0, 1, 0)) %>% 
  select(-GENE_EXP_at_least_one_match) 

ot_one <- coloc_results %>% 
  select(COLOCS_gene, VAR, twas_trait, Classif_Open_Targets_G2V_RS_value, 
         Open_targets_G2V_RS_value_gene) %>% 
  filter(!is.na(Open_targets_G2V_RS_value_gene)) %>% 
  group_by(VAR, twas_trait) %>% 
  summarise(open_targets_at_least_one_match = sum(COLOCS_gene == Open_targets_G2V_RS_value_gene)) %>% 
  ungroup() %>% 
  count(open_targets_at_least_one_match > 0) %>% 
  mutate(pct = n/sum(n)) %>% 
  slice(2)

ot_any_pairs <- coloc_results %>%
  select(COLOCS_gene, VAR, twas_trait, Classif_Open_Targets_G2V_RS_value, 
         Open_targets_G2V_RS_value_gene) %>% 
  # filter(!is.na(Open_targets_G2V_RS_value_gene)) %>% 
  group_by(VAR, twas_trait) %>% 
  summarise(open_targets_at_least_one_match = sum(COLOCS_gene == Open_targets_G2V_RS_value_gene)) %>% 
  ungroup() %>% 
  mutate(ot_any = if_else(open_targets_at_least_one_match > 0, 1, 0)) %>% 
  select(-open_targets_at_least_one_match)


ot_raw <- coloc_results %>% 
  select(COLOCS_gene, VAR, twas_trait, Classif_Open_Targets_G2V_RS_value, 
         Open_targets_G2V_RS_value_gene, 
         G2V_RS_value) 

ot_max = ot_raw %>% 
  group_by(VAR, twas_trait) %>% 
  filter(G2V_RS_value == max(G2V_RS_value)) %>% 
  select(VAR, twas_trait, ot_max_gene = Open_targets_G2V_RS_value_gene) %>% 
  distinct()

ot_max_table <- left_join(ot_raw, 
          ot_max, by = c("VAR", "twas_trait")) %>% 
  filter(!is.na(ot_max_gene)) %>% 
  group_by(VAR, twas_trait) %>% 
  arrange(COLOCS_gene) %>% 
  summarise(at_least_one_COLOCS_gene_is_ot_max = sum(COLOCS_gene == ot_max_gene)) %>% 
  ungroup() %>% 
  count(at_least_one_COLOCS_gene_is_ot_max > 0) %>% 
  mutate(pct = n/sum(n)) %>% 
  slice(2)

ot_max_pairs <- left_join(ot_raw, 
                          ot_max, by = c("VAR", "twas_trait")) %>% 
  group_by(VAR, twas_trait) %>% 
  arrange(COLOCS_gene) %>% 
  summarise(at_least_one_COLOCS_gene_is_ot_max = sum(COLOCS_gene == ot_max_gene)) %>% 
  ungroup() %>% 
  mutate(ot_max = if_else(at_least_one_COLOCS_gene_is_ot_max > 0, 1, 0)) %>% 
  select(-at_least_one_COLOCS_gene_is_ot_max)


set_df <- left_join(all_pairs, ot_max_pairs) %>% 
  left_join(ot_any_pairs) %>% 
  left_join(csq_pairs) %>% 
  left_join(gene_exp_pairs) %>% 
  mutate(id = str_c(VAR, "-",twas_trait)) %>% 
  select(-c(VAR, twas_trait))  %>% 
  filter(!if_all(where(is.numeric), ~is.na(.))) %>% 
  mutate(avail_ot_max = if_else(is.na(ot_max), 0, 1), 
         avail_csq = if_else(is.na(csq), 0, 1), 
         avail_gexp = if_else(is.na(gexp), 0, 1))

library(UpSetR)

set_df %>% 
  mutate(across(where(is.numeric), ~if_else(is.na(.x), 0, .x))) %>% 
  as.data.frame() %>% 
  upset(sets = c("ot_max", "ot_any", "csq", "gexp"), 
        order.by = 'freq')

set_df %>% 
  as.data.frame() %>% 
  upset(sets = c("avail_ot_max", "avail_csq", "avail_gexp"), 
        order.by = 'freq')


set_df %>% 
  mutate(across(where(is.numeric), ~if_else(is.na(.x), 0, .x))) %>% 
  as.data.frame() %>% 
  upset(sets = c("ot_max", "ot_any", "csq", "gexp"), 
        order.by = 'freq')

#restrict to genes in all 3 datasets
set_df %>% 
  filter(if_all(starts_with("avail"), ~ . == 1))%>% 
  as.data.frame() %>% 
  upset(sets = c("ot_max", "ot_any", "csq", "gexp"), 
        order.by = 'freq')

#restrict to genes in all 3 datasets
set_df %>% 
  filter(avail_ot_max == 1) %>% 
  mutate(across(where(is.numeric), ~if_else(is.na(.x), 0, .x))) %>% 
  as.data.frame() %>% 
  upset(sets = c("ot_max", "ot_any", "csq", "gexp"), 
        order.by = 'freq')


coloc_pct_table <- bind_rows(coloc_csq, 
          coloc_gene_exp, 
          ot_one, 
          ot_max_table) %>% 
  mutate(ref = c("csq", "gene_exp", 
                    "ot_any", "ot_max"), 
         method = "coloc") %>% 
  select(method, ref, n, pct)

write_tsv(coloc_pct_table, 
          "coloc_pct_table.tsv")

