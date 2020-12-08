library(tidyverse)
library(data.table)
library(janitor)
marg_results_credSet <- fread("../../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_FINEMAP_ref_file_TWASLoci.tsv") %>% 
  as_tibble()

gwas_to_twas <- fread("../../../data/BCX_GWAS_2020_known_variants/GWAS_sentinels_to_TWAS_genes.tsv") %>% 
  clean_names()

gwas_to_twas %>% distinct(twas_trait, unique_variant_id_b37)

#Construct level 1 gene list
data_level1 <- marg_results_credSet %>% 
  filter(!is.na(prob)) %>% 
  group_by(gene_name) %>% 
  select(chr:gene_name, model_r2, locus_name, phenotype, log10_regenie_p, marginal_significant) %>% 
  mutate(pheno_list = str_c(phenotype, collapse = ","), 
         n = n()) 

# 5049 distinct genes
level1_genes <- data_level1 %>% select(twas_gene = gene_name, 
                                          twas_trait = phenotype, model_r2)
#4888 distinct twas genes overlapping with sentinels. 
data_level1_GWAS <- inner_join(gwas_to_twas, level1_genes, 
                                  by = c("twas_trait", "twas_gene"))

data_level1_GWAS %>% distinct(twas_trait, unique_variant_id_b37)
###################
# Level 2 results #
###################

data_level2 <- marg_results_credSet %>% 
  filter(sig_high_PPCS == 1)

level2_genes <- data_level2 %>% select(twas_gene = gene_name, 
                          twas_trait = phenotype, 
                          dgn_chunk, model_r2)
#2,649 distinct TWAS genes
data_level2_GWAS <- inner_join(gwas_to_twas, level2_genes, 
           by = c("twas_trait", "twas_gene")) %>% 
  group_by(unique_variant_id_b37, twas_trait) %>% 
  mutate(n_sig_high_PPCS_genes = n()) %>% 
  ungroup()

data_level2_GWAS  %>% distinct(twas_trait, unique_variant_id_b37)

#Write distinct level 2 gene-variant pairs to compute correlations. 
# Also need DGN chunk to pull the correct expression script. 
# data_level2_GWAS %>%
#   distinct(unique_variant_id_b37, twas_gene, dgn_chunk, chr_gr_ch37) %>%
#   write_tsv("../../data/BCX_GWAS_2020_known_variants/twas_level2_varGene_pairs")

#Read in correlation results 
cor_files <- list.files("../../../data/BCX_GWAS_2020_known_variants/", 
                        pattern = "twas_level2_varGene_pairs_Correlation",
                        full.names = T)

l <- vector("list", length=length(cor_files))

for(i in 1:length(cor_files)){
  l[[i]] <- fread(cor_files[[i]]) %>% 
    as_tibble()
  
}

cor_data <- bind_rows(l) %>% 
  filter(!is.na(cor)) %>% 
  select(unique_variant_id_b37, twas_gene, cor)

data_level2_GWAS_cor <- data_level2_GWAS %>%  left_join(cor_data, 
                                                           by = c("unique_variant_id_b37",
                                                                  "twas_gene")) %>% 
  mutate(r2 = cor^2) %>% 
  filter(!is.na(r2))

#1,103 distinct genes
data_level2_GWAS_cor %>% 
  filter(r2 > .5) %>% 
  distinct(twas_gene)

#2,590 GWAS variants mapped to a gene
data_level2_GWAS_cor %>% 
  filter(r2 > .5) %>% 
  distinct(twas_trait, unique_variant_id_b37) 

data_level2_GWAS_cor %>% 
  filter(r2 > .5) %>% 
  count(twas_trait, unique_variant_id_b37) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram()

# 358 distinct genes to resolve with a conditional test. 
data_level2_GWAS_cor %>% 
  filter(r2 > .5) %>% 
  group_by(unique_variant_id_b37, twas_trait) %>% 
  mutate(n_genes = n()) %>% 
  filter(n_genes > 1) %>% 
  ungroup() %>% 
  distinct(twas_trait, twas_gene, chr_gr_ch37, dgn_chunk)

# %>% 
#   write_tsv("../../data/GWAS_conditional_on_TWAS/GWAS_on_TWAS_corGenes_GeneRef.tsv")


####################
# EDA section ######
####################
 
data_level2_GWAS_cor %>% 
  ggplot(aes(x = model_r2, y = r2)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  labs(caption = "all genes are significant")

#The inflation around 0 is somewhat reassuring, because it would make
# sense that the causal genes should be around some causal GWAS signal.
data_level2_GWAS_cor %>% 
  mutate(dis = abs(bp_gr_ch38 - gene_start))%>% 
  mutate(twas_sentinel = if_else(twas_gene == locus_name, 1, 0)) %>% 
  filter(dis < 1e6) %>% 
  ggplot(aes(x = dis, y = r2, color = as.factor(twas_sentinel))) + 
  geom_point(alpha = 0.3) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", 
             color = "red") + theme_bw() + 
  labs(x = "Distance (bp): 0-1MB",
       y = "R2 between dosage and expression") + 
  geom_abline(intercept = 0.8, 
              slope = (-.8/1e6), 
              linetype = "dashed",
              color = "blue")

data_level2_GWAS_cor %>% 
  group_by(unique_variant_id_b37, twas_trait) %>% 
  summarise(max_r2 = max(r2)) %>% 
  # filter(max_r2 > .15) %>% 
  ggplot(aes(x = max_r2)) + 
  geom_density() + 
  geom_vline(xintercept = .5)

#############
# Egene comp #
##############

#how does our variant-gene correlation measure identify egenes in our 
# data? 

egene <- fread("../../../data/BCX_eQTL_coloc_results/2019_07_09_ld_signif_eqtl_coloc_results_CLEAN.tsv") %>% 
  as_tibble()
egenes <- egene %>% 
  distinct(gene_name, UKB_TWAS_trait) %>% 
  mutate(is_egene = 1) %>% 
  rename(twas_gene = gene_name, 
         twas_trait = UKB_TWAS_trait)

egene %>% 
  distinct(condsig_var, trait)


egene_traits <- egenes %>% distinct(twas_trait)

twas_gene_list <- data_level2_GWAS_cor %>% 
  inner_join(egene_traits) %>% 
  filter(r2 > .5) %>% 
  #already distinct here. 
  distinct(twas_trait, unique_variant_id_b37, twas_gene) 


egene %>% 
  distinct(condsig_var, UKB_TWAS_trait, gene_name) %>% 
  group_by(condsig_var, UKB_TWAS_trait) %>% 
  summarise(eGenes = str_c(gene_name, collapse = ",")) %>% 
  rename(unique_variant_id_b37 = condsig_var, 
         twas_trait = UKB_TWAS_trait) %>% 
  left_join(twas_gene_list, .) %>% 
  # filter(!is.na(eGenes)) %>% 
  rowwise() %>% 
  mutate(check = if_else(str_detect(eGenes, twas_gene), 1, 0))
  

library(janitor)
tab <- data_level2_GWAS_cor %>% 
  inner_join(egene_traits, by = c("twas_trait")) %>% 
  distinct(twas_trait, twas_gene, r2) %>% 
  left_join(egenes) %>% 
  mutate(is_egene = if_else(is.na(is_egene), 0, is_egene), 
         r2_cutoff = if_else(r2 > 0.5, 1, 0)) %>% 
  tabyl(r2_cutoff, is_egene) 
#2 way table
tab
#chisq test
tab %>% chisq.test()

#row pct
tab %>% 
  adorn_totals(c("row", "col"))

#column pct
tab %>% 
  adorn_totals(c("row", "col")) %>% 
  adorn_percentages(c("col"))

###########
# GWAS gene #
#############

gwas_tab <- data_level2_GWAS_cor %>% 
  mutate(matches = map2_dbl(gene_symbol_s_for_most_serious_consequence, 
                            twas_gene, 
                            ~if_else(str_detect(.x, .y), 1, 0)), 
         r2_cutoff = if_else(r2 > .5, 1, 0)) %>% 
  tabyl(r2_cutoff, matches) 
gwas_tab%>% 
  adorn_totals(c("row", "col")) 

gwas_tab %>% chisq.test()

