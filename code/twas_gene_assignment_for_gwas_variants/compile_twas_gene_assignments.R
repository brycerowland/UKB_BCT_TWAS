library(tidyverse)
library(data.table)
library(janitor)
marg_results_credSet <- read_tsv("../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_FINEMAP_ref_file_TWASLoci.tsv") %>% 
  as_tibble()

gwas_to_twas <- read_tsv("../../data/BCX_GWAS_2020_known_variants/GWAS_sentinels_to_TWAS_genes.tsv") %>% 
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
cor_files <- list.files("../../data/BCX_GWAS_2020_known_variants/", 
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
  geom_vline(xintercept = .2)

twas_gene_assignment <- data_level2_GWAS_cor %>% filter(r2 > .2)

twas_gene_assignment %>% 
  distinct(associated_blood_index, associated_blood_index_class, 
         unique_variant_id_b37, 
         twas_trait,
         twas_gene, r2) %>% 
  write_tsv("gwas_variants_twas_gene_assignments.tsv")

#############
# Egene comp #
##############

#how does our variant-gene correlation measure identify egenes in our 
# data? 

egene <- fread("../../data/BCX_eQTL_coloc_results/2019_07_09_ld_signif_eqtl_coloc_results_CLEAN.tsv") %>% 
  as_tibble()

egene %>%
  distinct(condsig_var, gene_name, UKB_TWAS_trait) %>% 
  write_tsv("gwas_variants_coloc_gene_assignments.tsv")


egene_triples <- egene  %>% distinct(condsig_var, gene_name, UKB_TWAS_trait) %>% 
  mutate(assoc_id = str_c(condsig_var, "_", UKB_TWAS_trait))
egene_traits <- egene_triples %>% distinct(UKB_TWAS_trait) %>% 
  pull()

egene_triples %>% 
  distinct(condsig_var)

twas_triples <- twas_gene_assignment %>% 
  distinct(unique_variant_id_b37, twas_trait, twas_gene, 
           r2) %>% 
  filter(twas_trait %in% egene_traits) %>% 
  mutate(assoc_id = str_c(unique_variant_id_b37, "_", twas_trait))

library(ggVennDiagram)
venn_data <- list(
  coloc_vars = unique(egene_triples$assoc_id),
  twas_vars = unique(twas_triples$assoc_id)
)
saveRDS(venn_data, "../../data/finemapping/venn_data.rds")
egene_assoc_ids <- egene_triples %>%  distinct(assoc_id)
twas_assoc_ids <- twas_triples %>% distinct(assoc_id)

#Assoc_id assigned by both - coloc
inner_join(egene_triples, twas_assoc_ids, by = "assoc_id") %>% 
  write_tsv("../../data/finemapping/twas_v_coloc_ID_in_both_egene_triples.tsv")

#Assoc_id assigned by both - TWAS
inner_join(twas_triples, egene_assoc_ids, by = "assoc_id")  %>% 
  write_tsv("../../data/finemapping/twas_v_coloc_ID_in_both_twas_triples.tsv")


#Assoc_id assigned by only coloc
anti_join(egene_triples, twas_assoc_ids) %>% 
  write_tsv("../../data/finemapping/twas_v_coloc_ID_only_egene_twas_triples.tsv")

#Assoc_id assigned by only TWAS
anti_join(twas_triples, egene_assoc_ids)  %>% 
  write_tsv("../../data/finemapping/twas_v_coloc_ID_only_TWAS_twas_triples.tsv")

#How many associations are mapped using coloc vs. TWAS
ggVennDiagram(venn_data)

#Of the variants that are mapped in both, how often do the genes agree?
twas_list <- twas_triples %>% distinct(assoc_id, twas_gene) 

coloc_list <- egene_triples %>% distinct(assoc_id, gene_name)

#269 variants overlap
inner_join(twas_list,coloc_list) %>% distinct(assoc_id)

#167 variants have at least 1 gene agree. 219/269 = 80%
inner_join(twas_list, 
           coloc_list, 
           by = c("assoc_id", "twas_gene"="gene_name")) %>% 
  distinct(assoc_id)

## Of the variants that are identified by coloc and not by TWAS
# what is going on for the TWAS end? 

twas_var_list <- twas_list %>% 
  pull(assoc_id)

#158 associations mapped to 168 genes. 
twas_missed_coloc_mapped_list <- egene_triples %>% 
  filter(!(assoc_id %in% twas_var_list))


# egene %>% distinct(condsig_var, gene_name, UKB_TWAS_trait, 
#                    eqtl_sentinel, 
#                    eqtl_beta) %>% 
#   rename(condsig_var_b37 = condsig_var, 
#          eqtl_sentinel_b37 = eqtl_sentinel) %>% 
#   write_tsv("../../data/finemapping/coloc_assigned_variants_TWAS_misses.tsv")

  
twas_genes <- marg_results_credSet %>% 
  distinct(gene_name) %>% 
  pull(gene_name)

twas_vars <- gwas_to_twas %>% 
  distinct(unique_variant_id_b37) %>% 
  pull(unique_variant_id_b37)

#89 triplets aren't considered (model r2 < 0.05, not in DGN, 
# or no nearby TWAS genes with good models (all loci within +/- 1MB)).
# 87 associations
twas_missed_coloc_mapped_list %>% 
  filter(!((gene_name %in% twas_genes) & (condsig_var %in% twas_vars))) %>% 
  distinct(assoc_id)


#Specifically, 36 vars not mapped to a TWAS gene. 13 don't have genes with DGN expression, 
# 23 have poorly trained TWAS models.
twas_missed_coloc_mapped_list %>% filter(!(gene_name %in% twas_genes)) %>% distinct(assoc_id)

#and, 52 not mapped to nearby TWAS genes. 
twas_missed_coloc_mapped_list %>% filter((gene_name %in% twas_genes), 
                                         !(condsig_var %in% twas_vars)) %>% distinct(assoc_id)

#74 associations remain.
well_predicted_missed_genes <- twas_missed_coloc_mapped_list %>% 
  filter(gene_name %in% twas_genes, 
         condsig_var %in% twas_vars)
well_predicted_missed_genes %>% distinct(assoc_id) %>% nrow()

#Of 74 remaining associations, 42 are to genes not TWAS significant. 
marg_results_credSet %>% 
  inner_join(well_predicted_missed_genes, 
             by = c("gene_name", "phenotype"="UKB_TWAS_trait")) %>% 
  filter(marginal_significant == 0) %>% 
  distinct(assoc_id)

# and 9 are not included in the FINEMAP credible set.
marg_results_credSet %>% 
  inner_join(well_predicted_missed_genes, 
             by = c("gene_name", "phenotype"="UKB_TWAS_trait")) %>% 
  filter((marginal_significant == 1 & is.na(prob))) %>% 
  distinct(assoc_id)
 
# 1 - 138/158 variants = 13% missed by TWAS correlation. Other filters
# feel very sensible to me. 


sig_well_predicted_twas_misses <- marg_results_credSet %>% 
  inner_join(well_predicted_missed_genes, 
             by = c("gene_name", "phenotype"="UKB_TWAS_trait")) %>% 
 filter(marginal_significant == 1, !is.na(prob))%>% 
  filter(prob > 0.5) %>% 
  select(unique_variant_id_b37 = condsig_var, 
         twas_trait = phenotype, 
         twas_gene = gene_name, prob)

sig_well_predicted_twas_misses %>% distinct(unique_variant_id_b37, twas_trait)


#Remaining 70 genes do not pass the correlation
# filter. 58 variants. 
sig_well_predicted_twas_misses %>% 
  distinct(unique_variant_id_b37)

data_level2_GWAS_cor %>% 
  inner_join(sig_well_predicted_twas_misses, 
             by = c("unique_variant_id_b37", 
                    "twas_trait", 
                    "twas_gene")) %>% 
  ggplot(aes(x = r2)) + geom_histogram()
data_level2_GWAS_cor %>% 
  inner_join(well_predicted_missed_genes)

#################
# What about this comparison with pheno class
# instead of trait? 
###################

egene_triples_class <- egene %>%
  mutate(pheno_class = case_when(
    str_detect(UKB_TWAS_trait,"baso|eos|lymph|neutr|mono|white") ~ "WBC",
    str_detect(UKB_TWAS_trait, "plate") ~ "PLT",
    TRUE ~ "RBC"
  )) %>% 
  distinct(condsig_var, gene_name, UKB_TWAS_trait, pheno_class) %>% 
  mutate(assoc_id = str_c(condsig_var, "_", pheno_class))


twas_triples_class <- twas_gene_assignment %>% 
  distinct(unique_variant_id_b37, twas_trait, twas_gene, associated_blood_index_class) %>% 
  filter(twas_trait %in% egene_traits) %>% 
  mutate(assoc_id = str_c(unique_variant_id_b37, "_", associated_blood_index_class))


venn_data_class <- list(
  coloc_vars = unique(egene_triples_class$assoc_id),
  twas_vars = unique(twas_triples_class$assoc_id)
)

#How many associations are mapped using coloc vs. TWAS
ggVennDiagram(venn_data_class)

#Of the variants that are mapped in both, how often do the genes agree?
twas_list <- twas_triples_class %>% 
  distinct(assoc_id, twas_gene) 

coloc_list <- egene_triples_class %>% 
  distinct(assoc_id, gene_name)

#177 variants overlap
inner_join(twas_list, 
           coloc_list) %>% distinct(assoc_id)

#167 variants have at least 1 gene agree. 145/177 = 80%
inner_join(twas_list, 
           coloc_list, 
           by = c("assoc_id", "twas_gene"="gene_name")) %>% 
  distinct(assoc_id)

## Of the variants that are identified by coloc and not by TWAS
# what is going on for the TWAS end? 
twas_var_list <- twas_list %>% 
  pull(assoc_id)

#174 variants mapped to 226 genes. 
twas_missed_coloc_mapped_list <- egene_triples_class %>% 
  filter(!(assoc_id %in% twas_var_list))
twas_missed_coloc_mapped_list %>% distinct(assoc_id) %>% nrow()

twas_genes <- marg_results_credSet %>% 
  distinct(gene_name) %>% 
  pull(gene_name)

twas_vars <- gwas_to_twas %>% 
  distinct(unique_variant_id_b37) %>% 
  pull(unique_variant_id_b37)

#97 triplets aren't considered (model r2 < 0.05, not in DGN, 
# or no nearby TWAS genes with good models (all loci within +/- 1MB)).
twas_missed_coloc_mapped_list %>% 
  filter(!((gene_name %in% twas_genes) & (condsig_var %in% twas_vars)))

#129 remain.
well_predicted_missed_genes <- twas_missed_coloc_mapped_list %>% 
  filter(gene_name %in% twas_genes, 
         condsig_var %in% twas_vars)

#Of 129 remaining genes. 48 are not marginal TWAS significant. 
marg_results_credSet %>% 
  inner_join(well_predicted_missed_genes, 
             by = c("gene_name", "phenotype"="UKB_TWAS_trait")) %>% 
  filter(marginal_significant == 0)

# and 12 are not included in the FINEMAP credible set.
marg_results_credSet %>% 
  inner_join(well_predicted_missed_genes, 
             by = c("gene_name", "phenotype"="UKB_TWAS_trait")) %>% 
  filter(marginal_significant == 1, is.na(prob))

# So 70% of genes have decent reasons to not be mapped. 
# 56 variants remain.
sig_well_predicted_twas_misses <- marg_results_credSet %>% 
  inner_join(well_predicted_missed_genes, 
             by = c("gene_name", "phenotype"="UKB_TWAS_trait")) %>% 
  filter(marginal_significant == 1, !is.na(prob)) %>% 
  select(unique_variant_id_b37 = condsig_var, 
         twas_trait = phenotype, 
         twas_gene = gene_name, prob) 



#Remaining 69 genes do not pass the correlation
# filter.
data_level2_GWAS_cor %>% 
  inner_join(sig_well_predicted_twas_misses, 
             by = c("unique_variant_id_b37", 
                    "twas_trait", 
                    "twas_gene")) %>% 
  ggplot(aes(x = r2)) + geom_histogram()






####################

twas_egene_merge <- egene %>% 
  distinct(condsig_var, UKB_TWAS_trait, gene_name) %>% 
  group_by(condsig_var, UKB_TWAS_trait) %>% 
  summarise(eGenes = str_c(gene_name, collapse = ",")) %>% 
  rename(unique_variant_id_b37 = condsig_var, 
         twas_trait = UKB_TWAS_trait) %>% 
  left_join(twas_gene_list, .) 

twas_egene_merge %>% 
  filter(!is.na(eGenes)) %>%
  rowwise() %>% 
  mutate(check = if_else(str_detect(eGenes, twas_gene), 1, 0)) %>% 
  ungroup() %>% 
  summarise(prop_match = mean(check))
  

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

