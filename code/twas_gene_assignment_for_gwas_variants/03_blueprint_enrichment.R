library(tidyverse)
library(kableExtra)

source("../functional_annotation_TWAS_results/twas_trait_to_cell_type.R")

coloc_v2g <- read_tsv("gwas_variants_coloc_gene_assignments.tsv")
coloc_phenos <- coloc_v2g %>% 
  pull(UKB_TWAS_trait) %>% 
  unique()

twas_v2g <- read_tsv("gwas_variants_twas_gene_assignments.tsv") 
# %>% 
#   filter(twas_trait %in% coloc_phenos)

both_coloc <- read_tsv("../../data/finemapping/twas_v_coloc_ID_in_both_egene_triples.tsv")
both_twas <- read_tsv("../../data/finemapping/twas_v_coloc_ID_in_both_twas_triples.tsv")
only_twas <- read_tsv("../../data/finemapping/twas_v_coloc_ID_only_TWAS_twas_triples.tsv")
only_coloc <- read_tsv("../../data/finemapping/twas_v_coloc_ID_only_egene_twas_triples.tsv")

ensg_id_key <- read_tsv("/proj/yunligrp/users/bryce/twas/reference_files/gencode_build38_genes.txt", 
                        col_names = c("chr", "start_pos", "end_pos", 
                                      "ensg_id", "gene_name")) %>% 
  rowwise() %>% 
  mutate(ensg_id = str_remove(ensg_id, "\\.[0-9]*$")) %>% 
  ungroup() %>% 
  select(ensg_id, gene_name)


cell_type_specific_files <- list.files("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon", 
                                       pattern = "CellTypeSpecificGene", 
                                       full.names = T)
l <- vector("list", length = length(cell_type_specific_files))

for(i in 1:length(cell_type_specific_files)){
  cell_type <- cell_type_specific_files[[i]] %>% 
    basename() %>% 
    str_remove("CellTypeSpecificGene_") %>% 
    str_remove(".txt")
  
  l[[i]] <- read_tsv(cell_type_specific_files[[i]]) %>% 
    mutate(cell_type = cell_type)
}


cts_data <- bind_rows(l) 

cts_check <- function(.ensg_id, .cell_type){
  cts_data %>% 
    filter(ENSEMBL_GENEID == .ensg_id) %>% 
    distinct(cell_type) %>% 
    pull(cell_type)
}



#Both TWAS
twas_cts <- left_join(twas_v2g, ensg_id_key, 
          by = c("twas_gene" = "gene_name")) %>% 
  rowwise() %>% 
  mutate(twas_cell_type = list(twas_trait_to_cell_type(twas_trait)), 
         cts = list(cts_check(ensg_id, twas_cell_type)), 
         cts_str = str_c(cts, collapse = ","), 
         cts_exact_match = if_else(any(twas_cell_type %in% cts), 1, 0), 
         cts_some_match = if_else(cts_str != "", 1, 0)) %>% 
  ungroup() 

twas_cts_var_results <- twas_cts %>% 
  select(unique_variant_id_b37, twas_trait, cts_exact_match, cts_some_match) %>% 
  group_by(unique_variant_id_b37, twas_trait) %>% 
  summarise(s_exact = sum(cts_exact_match), 
            s_some = sum(cts_some_match)) %>% 
  ungroup() %>% 
  summarise(n_var_exact_match = sum(s_exact > 0), 
            pct_var_exact_match = mean(s_exact > 0), 
            n_var_some_match = sum(s_some > 0), 
            pct_var_some_match = mean(s_some > 0))

twas_formatted <- twas_cts_var_results %>% 
  pivot_longer(cols = everything(),
               names_to = c("measure", "ref"),
               names_pattern = "(.*)_var_(.*)") %>% 
  pivot_wider(id_cols = ref, 
              names_from = "measure", 
              values_from = "value") %>% 
  mutate(method = "twas", 
         ref = str_c("blueprint_", ref))


#Both coloc
coloc_cts <- left_join(coloc_v2g, ensg_id_key, 
                      by = c("gene_name")) %>% 
  rowwise() %>% 
  mutate(twas_cell_type = twas_trait_to_cell_type(UKB_TWAS_trait), 
         cts = list(cts_check(ensg_id, twas_cell_type)), 
         cts_str = str_c(cts, collapse = ","), 
         cts_exact_match = if_else(twas_cell_type %in% cts, 1, 0), 
         cts_some_match = if_else(cts_str != "", 1, 0)) %>% 
  ungroup()

coloc_cts_var_results <- coloc_cts %>% 
  select(condsig_var, UKB_TWAS_trait, cts_exact_match, cts_some_match) %>% 
  group_by(condsig_var, UKB_TWAS_trait) %>% 
  summarise(s_exact = sum(cts_exact_match), 
            s_some = sum(cts_some_match)) %>% 
  ungroup() %>% 
  summarise(n_var_exact_match = sum(s_exact > 0), 
            pct_var_exact_match = mean(s_exact > 0), 
            n_var_some_match = sum(s_some > 0), 
            pct_var_some_match = mean(s_some > 0))

coloc_formatted <- coloc_cts_var_results %>% 
  pivot_longer(cols = everything(),
               names_to = c("measure", "ref"),
               names_pattern = "(.*)_var_(.*)") %>% 
  pivot_wider(id_cols = ref, 
              names_from = "measure", 
              values_from = "value") %>% 
  mutate(method = "coloc", 
         ref = str_c("blueprint_", ref))

bind_rows(twas_formatted, 
          coloc_formatted) %>% 
  write_tsv("soranzo_group_results/pct_blueprint_table_all_phenos.tsv")

# Only TWAS
only_twas_cts <- left_join(only_twas, ensg_id_key, 
                           by = c("twas_gene" = "gene_name")) %>% 
  rowwise() %>% 
  mutate(twas_cell_type = twas_trait_to_cell_type(twas_trait), 
         cts = list(cts_check(ensg_id, twas_cell_type)), 
         cts_str = str_c(cts, collapse = ","), 
         cts_exact_match = if_else(twas_cell_type %in% cts, 1, 0), 
         cts_some_match = if_else(cts_str != "", 1, 0)) %>% 
  ungroup()

only_twas_cts_var_results <- only_twas_cts %>% 
  select(assoc_id, cts_exact_match, cts_some_match) %>% 
  group_by(assoc_id) %>% 
  summarise(s_exact = sum(cts_exact_match), 
            s_some = sum(cts_some_match)) %>% 
  summarise(n_var_exact_match = sum(s_exact > 0), 
            pct_var_exact_match = mean(s_exact > 0), 
            n_var_some_match = sum(s_some > 0), 
            pct_var_some_match = mean(s_some > 0))


only_twas_cts %>% 
  summarise(n_exact_match = sum(cts_exact_match), 
            pct_exact_match = mean(cts_exact_match), 
            n_some_match = sum(cts_some_match), 
            pct_some_match = mean(cts_some_match))

# Only coloc
only_coloc_cts <- left_join(only_coloc, ensg_id_key, 
                           by = c("gene_name" = "gene_name")) %>% 
  rowwise() %>% 
  mutate(twas_cell_type = twas_trait_to_cell_type(UKB_TWAS_trait), 
         cts = list(cts_check(ensg_id, twas_cell_type)), 
         cts_str = str_c(cts, collapse = ","), 
         cts_exact_match = if_else(twas_cell_type %in% cts, 1, 0), 
         cts_some_match = if_else(cts_str != "", 1, 0)) %>% 
  ungroup()

only_coloc_cts_var_results <- only_coloc_cts %>% 
  select(assoc_id, cts_exact_match, cts_some_match) %>% 
  group_by(assoc_id) %>% 
  summarise(s_exact = sum(cts_exact_match), 
            s_some = sum(cts_some_match)) %>% 
  summarise(n_var_exact_match = sum(s_exact > 0), 
            pct_var_exact_match = mean(s_exact > 0), 
            n_var_some_match = sum(s_some > 0), 
            pct_var_some_match = mean(s_some > 0))


only_coloc_cts %>% 
  summarise(n_exact_match = sum(cts_exact_match), 
            pct_exact_match = mean(cts_exact_match), 
            n_some_match = sum(cts_some_match), 
            pct_some_match = mean(cts_some_match))
library(kableExtra)

#Assemble results
bind_rows(both_twas_cts_var_results, both_coloc_cts_var_results,
          only_twas_cts_var_results, only_coloc_cts_var_results) %>% 
  mutate(variant_set = c("Both TWAS", "Both Coloc",
                         "ONLY TWAS", "Only Coloc")) %>% 
  relocate(variant_set) %>% 
  rename(`Association Set` = variant_set, 
         `Exact Match (n)` = n_var_exact_match,
         `Exact Match (%)` = pct_var_exact_match, 
         `Some Match (n)` = n_var_some_match, 
         `Some Match (%)` = pct_var_some_match) %>% 
  kable("latex", booktabs = T, digits = 3)


%>% 
  View()

