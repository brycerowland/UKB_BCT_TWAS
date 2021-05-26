library(tidyverse)
library(LDlinkR)
setwd("/proj/yunligrp/users/bryce/twas/UKB/code/twas_gene_assignment_for_gwas_variants")

#Read in missed GWAS variant file which contains the eQTL information identified by 
# coloc. Note that some GWAS variants (most) have more than one eQTL. 
#
#Also, the missed_vars file are acutally just all of the coloc results. 

missed_vars <- read_tsv("../../data/finemapping/coloc_assigned_variants_TWAS_misses.tsv")%>% 
  rowwise() %>% 
  mutate(chr = str_split(condsig_var_b37, ":")[[1]][1]) %>% 
  ungroup()


# Read in both-twas and both-coloc assignments 
both_coloc <- read_tsv("../../data/finemapping/twas_v_coloc_ID_in_both_egene_triples.tsv") %>% 
  rename(unique_variant_id_b37 = condsig_var, 
         twas_trait = UKB_TWAS_trait)
both_twas <- read_tsv("../../data/finemapping/twas_v_coloc_ID_in_both_twas_triples.tsv") %>% 
  rename(gene_name = twas_gene)

#167 associations completely agree between methods. Therefore, for these we have a 
# gene expression prediction model that is TWAS significant and correlated with 
# GWAS variant dosage, and an eQTL which is in high LD with the causal variant, all 
# for the same genes. 
#These correspond to 172 gene assignments. 
same_assignments <-inner_join(both_twas, 
           both_coloc, 
           by = c("unique_variant_id_b37",
                  "twas_trait",
                  "assoc_id",
                  "gene_name")) 

# The 172 gene assignments correspond to 197 eqtl's colocalizing with the GWAS variants. 
#
# Fortunately, I had already pulled their LD buddies since they're in the
# missed data set. 
egene_vars <-read_tsv("../../data/BCX_eQTL_coloc_results/2019_07_09_ld_signif_eqtl_coloc_results_CLEAN.tsv") %>% 
  select(twas_trait = UKB_TWAS_trait, 
         gene_name, 
         unique_variant_id_b37 = condsig_var, 
         eqtl_sentinel_b37 = eqtl_sentinel, 
         eqtl_beta) %>% 
  inner_join(same_assignments, by =c("twas_trait", "gene_name", "unique_variant_id_b37")) %>% 
  distinct() %>% 
  rowwise() %>%
  mutate(
    ldlink_var_id = str_split(eqtl_sentinel_b37, "_")[[1]][1] %>% 
      str_c("chr", .), 
    chr = str_split(eqtl_sentinel_b37, ":")[[1]][1]
  ) %>% 
  ungroup()

#For each variant-trait:gene pair, we want to know if there is at least 1 eqtl (or LD buddy of
# the eQTL in the gene expression prediction model of the assigned gene. 
ld_link_missed_vars <- missed_vars %>% 
  rowwise() %>%
  mutate(
    ldlink_var_id = str_split(eqtl_sentinel_b37, "_")[[1]][1] %>% 
      str_c("chr", .)
  ) %>% 
  ungroup()

eqtl_varid_list <- ld_link_missed_vars %>% 
  distinct(ldlink_var_id)


get_ld_buddies <- function(eqtl, r2_cutoff){
  file_name <- str_c("../../data/finemapping/coloc_assigned_variants_TWAS_model_liftover/eqtl_ld_buddies/",
                     eqtl, ".txt")
  
  if(file.exists(file_name)){
    read_tsv(file_name, 
             skip = 1, col_names = F) %>% 
      select(c(3, 8)) %>% 
      rename(ld_buddy_coords = X3, 
             R2 = X8) %>% 
      filter(R2 > r2_cutoff) %>% 
      pull(ld_buddy_coords) %>% return()
  }
  else{
    return(NULL)
  }
}

get_model_vars <- function(chr, gene_name){
  file_name <- str_c("../../data/finemapping/coloc_assigned_variants_TWAS_model_liftover/model_b37_ids_chr", chr, "_", gene_name)
  
  if(file.exists(file_name)){
    read_tsv(file_name, 
             col_names = "model_vars") %>% 
      rowwise() %>% 
      mutate(split = str_split(model_vars, ":"), 
             new_id = str_c("chr", split[1], ":",
                            split[2])) %>% 
      pull(new_id) %>% 
      return()
  }else{
    return(NULL)
  }
  
}


var_check_results <- ld_link_missed_vars %>% 
  rowwise() %>% 
  mutate(ld_buddies_list = list(get_ld_buddies(ldlink_var_id,
                                               0.5)), 
         model_vars_list = list(get_model_vars(chr, gene_name)),
         eqtl_in_model = if_else(ldlink_var_id %in% model_vars_list, 
                                 1, 0),
         n_buddies_in_model = sum(ld_buddies_list %in% model_vars_list),
         buddies_in_model = if_else(n_buddies_in_model > 0,
                                    1,0),
         either = if_else(eqtl_in_model == 1 | 
                            buddies_in_model == 1, 1, 0)) %>% 
  ungroup()

egene_var_check_results <- egene_vars %>% 
  rowwise() %>% 
  mutate(ld_buddies_list = list(get_ld_buddies(ldlink_var_id,
                                               0.5)), 
         model_vars_list = list(get_model_vars(chr, gene_name)),
         eqtl_in_model = if_else(ldlink_var_id %in% model_vars_list, 
                                 1, 0),
         n_buddies_in_model = sum(ld_buddies_list %in% model_vars_list),
         buddies_in_model = if_else(n_buddies_in_model > 0,
                                    1,0),
         either = if_else(eqtl_in_model == 1 | 
                            buddies_in_model == 1, 1, 0)) %>% 
  ungroup()


#Only need to run below calculation once. 

# setwd("/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/coloc_assigned_variants_TWAS_model_liftover/eqtl_ld_buddies")
# LDproxy_batch(eqtl_varid_list, 
#               pop = c("CEU", "TSI", "FIN", "GBR", "IBS"),
#               token = Sys.getenv("LDLINK_TOKEN"))


#Write a function that given an eqtl ID in LDLink format, read in the LD file, 
# filter it to a specific (R2 > x) threshold, and return LD buddies. 


var_check <- var_check_results %>% 
  mutate(assoc_id = str_c(condsig_var_b37, "_",
                          UKB_TWAS_trait)) 

either_score <- var_check %>% 
  group_by(assoc_id) %>% 
  summarise(s = sum(either))

either_score %>% 
  mutate(binary_s = if_else(s > 0 , 1, 0)) %>% 
  summarise(mean(binary_s))


#Subset check

subset_genes <- data_level2_GWAS_cor %>% 
  inner_join(sig_well_predicted_twas_misses, 
             by = c("unique_variant_id_b37", 
                    "twas_trait", 
                    "twas_gene")) %>% 
  select(unique_variant_id_b37, twas_gene, 
         twas_trait, r2)


subset_gene_var_check <- inner_join(var_check_results, subset_genes, 
           by = c("condsig_var_b37"="unique_variant_id_b37",
                  "gene_name"="twas_gene",
                  "UKB_TWAS_trait"="twas_trait")) %>% 
  mutate(assoc_id = str_c(condsig_var_b37, "_",
                          UKB_TWAS_trait)) 



either_score <- subset_gene_var_check %>% 
  group_by(assoc_id) %>% 
  summarise(s = sum(either))

prediction_results <- read_tsv("../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_ref_file_TWASLoci.tsv")%>% 
  select(gene_name, model_r2, log10_regenie_p, phenotype) %>% 
  distinct()
prediction_results 

plot_data <- inner_join(subset_gene_var_check, either_score) %>% 
  left_join(prediction_results, by = c("gene_name", 
                                       "UKB_TWAS_trait"="phenotype")) %>% 
  rowwise() %>% 
  mutate(n_model_vars = length(model_vars_list), 
         n_ld_buddies = length(ld_buddies_list)) %>% 
  distinct(assoc_id, gene_name, model_r2, s, n_model_vars, n_ld_buddies,
           log10_regenie_p) %>% 
  mutate(binary_s = if_else(s > 0, 1, 0)) %>% 
  ungroup()

plot_data %>% 
  distinct(assoc_id, gene_name, binary_s) %>% 
  group_by(assoc_id) %>% 
  summarise(sum = sum(binary_s)) %>% 
  summarise(mean(sum > 0))

 
plot_data%>% 
  ggplot(aes(x = as.factor(binary_s), y = model_r2)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point()

plot_data%>% 
  ggplot(aes(x = as.factor(binary_s), y = log10_regenie_p)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point()

plot_data%>% 
  ggplot(aes(x = as.factor(binary_s), y = n_model_vars)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point()

plot_data%>% 
  ggplot(aes(x = as.factor(binary_s), y = n_ld_buddies)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point() + 
  theme_bw() + 
  labs(x = "eQTL (or proxy) in TWAS model", 
       y = "# of LD buddies for eQTL")


comp_data <- bind_rows(egene_var_check_results   %>% 
  mutate(var_list = "Variant-Gene assignment\nAgree in Both"), 
subset_gene_var_check %>% 
  rename(unique_variant_id_b37 = condsig_var_b37, 
         twas_trait = UKB_TWAS_trait) %>% 
  mutate(var_list = "Variant missed by TWAS\nBut in Coloc")) %>% 
  left_join(prediction_results, by = c("gene_name", 
                                       "twas_trait"="phenotype"))%>% 
  rowwise() %>% 
  mutate(n_model_vars = length(model_vars_list), 
         n_ld_buddies = length(ld_buddies_list))

#Right off the bat, in 92% of mappings which agree between
# coloc and TWAS the eQTL is in the GEXp model, compared to 
# ~40% for those missed by TWAS. For the 60% where that's not true,
# we can assume that's why they're missed - they do not identify
# eQTL in LD with the GWAS variant we are trying to assign in their
# prediction models (however, they must identify some other signal)
# at the TWAS locus because the gene is TWAS significant. 
comp_data %>% 
  group_by(var_list) %>% 
  summarise(mean(either))

#What about genes with either == 1? Why are they missed? 

comp_data%>% 
  ggplot(aes(x = var_list, y = n_model_vars)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point()

#No real difference in eQTL effect size. 
comp_data %>% 
  filter(either == 1) %>%
  ggplot(aes(x = var_list, y = abs(eqtl_beta))) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point(alpha = 0.25)


#Definitely a difference. Genes missed by TWAS 
# are more poorly predicted in general. 
comp_data %>% 
  filter(either == 1) %>% 
  ggplot(aes(x = var_list, y = model_r2)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point(alpha = 0.25) 

#No real difference
comp_data %>% 
  filter(either == 1) %>% 
  ggplot(aes(x = var_list, y = log10_regenie_p)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point() + 
  coord_cartesian(ylim = c(0,100))

#Best evidence so far. 
comp_data %>% 
  mutate(prop_buddies_in_model = n_buddies_in_model/n_model_vars) %>% 
  # filter(either == 1) %>% 
  ggplot(aes(x = var_list, y = prop_buddies_in_model)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, .75)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  coord_flip()+ 
  labs(x = "Variant List", 
       y = "Proportion of TWAS model variants that are LD buddies
       for co-localized eQTL.") + 
  theme(text = element_text(size = 14))
  

comp_data %>% 
  filter(either == 1) %>% 
  ggplot(aes(x = n_ld_buddies, 
             y = n_model_vars)) + 
  geom_point(aes(color = var_list)) + 
  geom_abline()

comp_data %>% 
  ggplot(aes(x = n_ld_buddies, 
             y = model_r2)) + 
  geom_point(aes(color = var_list)) + 
  geom_abline()

library(janitor)
comp_data %>% 
  filter(either == 1) %>% 
  mutate(more_ld_buddies_than_model_vars = 
           if_else(n_ld_buddies > n_model_vars, 
                   "more_ld_buddies", "less"))  %>% 
  tabyl(var_list, 
        more_ld_buddies_than_model_vars) %>% 
  adorn_percentages("row")



