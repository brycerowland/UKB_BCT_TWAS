library(tidyverse)

files <- list.files("../../../data/conditional_analysis/REGENIE/analysis_results/", pattern = "_CALoci$", 
                    full.names = T)

l <- vector("list", length = length(files))

for(i in 1:length(files)){
  l[[i]] <- read_tsv(files[[i]], 
                     col_types = cols(
                       chr = col_double(),
                       start_pos = col_double(),
                       end_pos = col_double(),
                       gene_name = col_character(),
                       dgn_chunk = col_double(),
                       en_fit = col_double(),
                       model_r2 = col_double(),
                       cv_r2 = col_logical(),
                       log10_regenie_p = col_double(),
                       locus_name = col_character(),
                       phenotype = col_character(),
                       log10_conditional_p = col_double(),
                       conditional_locus = col_character()
                     ))
}




ca_results <- bind_rows(l) %>% 
  #Exclude HLA.
  filter(!(chr == 6 & 
             start_pos > 24477797 &
             end_pos < 37448354)) %>% 
  filter(log10_regenie_p > 6.763558) %>% 
  mutate(pheno_class = case_when(
    str_detect(phenotype,"baso|eos|lymph|neutr|mono|white") ~ "WBC", 
    str_detect(phenotype, "plate") ~ "PLT", 
    TRUE ~ "RBC"
  ))

ca_results %>% filter(log10_conditional_p > -1 *log10(0.05 / 11759)) %>% 
  count(gene_name) %>% 
  count(n > 1)


write_tsv(ca_results,
          "../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_CA_results.tsv")


gwas_vars <- read_tsv("../../../data/BCX_GWAS_2020_known_variants/Vuckovic_SuppTables_clean.tsv") %>% 
  distinct(unique_variant_id_b37, associated_blood_index, .keep_all = T)  %>% 
  mutate(unique_variant_id_alt_b38 = map_chr(unique_variant_id_b38,
                                             ~str_c(str_split(.x, "_")[[1]][c(1,3,2)],
                                                    collapse = "_")))
gwas_to_twas_pheno <- read_csv("/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/GWAS_to_TWAS.csv")

gwas <- left_join(gwas_vars, gwas_to_twas_pheno, by = c("associated_blood_index"="GWAS_trait"))

#For each gene start pos and end pos, we check to see if  there are any
# variants within the +/- 1MB window around the gene. 
#
# Returns the number of genes in the window. 
is_novel <- function(.chr, .start_pos, .end_pos){
  vars <- gwas %>% 
    filter(chr_gr_ch37 == .chr, 
           bp_gr_ch38 >  (.start_pos - 1e6), 
           bp_gr_ch38 < (.end_pos + 1e6)) %>% 
    nrow()
  return(vars)
}

pull_novel <- function(.chr, .start_pos, .end_pos){
  vars <- gwas %>% 
    filter(chr_gr_ch37 == .chr, 
           bp_gr_ch38 >  (.start_pos - 1e6), 
           bp_gr_ch38 < (.end_pos + 1e6)) %>% 
    select(unique_variant_id_alt_b38, var_reported_pheno = TWAS_trait, 
           bp_gr_ch38)
  return(vars)
}

# Same analysis as above, however we restrict the variants to be a 
# member of the BCT index class. 
is_novel_by_class <- function(.pheno_class, .chr, .start_pos, .end_pos){
  vars <- gwas %>% 
    filter(associated_blood_index_class == .pheno_class, 
           chr_gr_ch37 == .chr, 
           bp_gr_ch38 >  (.start_pos - 1e6), 
           bp_gr_ch38 < (.end_pos + 1e6)) %>% 
    nrow()
  return(vars)
}

ca_results_var_window <- ca_results %>% 
  mutate(vars_in_window = pmap_dbl(list(.chr = chr, 
                                        .start_pos = start_pos, 
                                        .end_pos = end_pos), is_novel)) 

ca_results_var_window %>% 
  count(vars_in_window)

ca_results_var_window_byClass <- ca_results %>% 
  mutate(vars_in_window = pmap_dbl(list(.chr = chr, 
                                        .start_pos = start_pos, 
                                        .end_pos = end_pos), is_novel)) %>% 
  mutate(vars_in_window_byClass = pmap_dbl(list(.pheno_class = pheno_class, 
                                        .chr = chr, 
                                        .start_pos = start_pos, 
                                        .end_pos = end_pos), is_novel_by_class)) 

sentinel_genes <- ca_results_var_window_byClass %>% 
  distinct(locus_name)

ca_results_var_window_byClass %>% 
  inner_join(sentinel_genes, by = c("gene_name"="locus_name")) %>% 
  select(gene_name, vars_in_window) %>% 
  group_by(gene_name) %>% 
  summarise(n = max(vars_in_window)) %>% 
  filter(n > 0)

%>% 
  distinct(locus_name, phenotype)

%>% 
  count(vars_in_window_byClass) %>% 
  print(n = 100)

# write_tsv(ca_results_var_window, 
#           path = "../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_novel_loci_info.tsv")
# 
# write_tsv(ca_results_var_window_byClass,
#           path = "../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_novel_lociByClass_info.tsv")

ca_results_var_window <- read_tsv("../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_novel_loci_info.tsv")
ca_results_var_window_byClass <- read_tsv("../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_novel_lociByClass_info.tsv")

### Variants for any Blood Cell Trait

#9 loci with no known variants in window of any genes at locus.
# (all are only marginal significant gene at locus)
novel_loci_genes <- ca_results_var_window %>% 
  group_by(conditional_locus, chr, phenotype, pheno_class) %>% 
  summarise(var_at_locus = sum(vars_in_window), genes_in_locus = n()) %>% 
  filter(var_at_locus == 0) %>% 
  ungroup() %>% 
  select(conditional_locus, phenotype)


inner_join(ca_results_var_window, novel_loci_genes, 
           by = c("conditional_locus", "phenotype")) %>% 
  write_tsv("../../../data/functional_annotation_TWAS_results/no_known_gwas_variants_genes.tsv")

get_marg_beta <- function(gene_name, pheno){
  read_delim(str_c("../../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_", 
             pheno, "_bonferroniSigGenes"),
             delim = " ") %>% 
    filter(GENE == gene_name) %>% 
    pull(BETA) %>% 
    return()
}
get_cond_beta <- function(gene_name, pheno){
  read_delim(str_c("../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_conditionalAnalysis_",
                   pheno, "_results"), 
             delim = " ") %>% 
    filter(ID == gene_name) %>% 
    pull(BETA) %>% 
    return()
}

library(kableExtra)
read_csv("copy_novel_genes.csv") %>% 
  select(gene_name, model_r2, log10_regenie_p, 
         log10_conditional_p, phenotype) %>% 
  rowwise() %>% 
  mutate(marg_beta = get_marg_beta(gene_name, 
                              phenotype), 
         cond_beta = get_cond_beta(gene_name, 
                                   phenotype)) %>% 
  ungroup() %>% 
  mutate(phenotype = phenotype %>% 
           str_replace_all("_", " ") %>% 
           str_to_title()) %>% 
  select(Gene = gene_name,
         `Exp. Model R2` = model_r2, 
         Phenotype = phenotype, 
         `TWAS Beta` = marg_beta, 
         `-log10 TWAS p` = log10_regenie_p) %>% 
  write_tsv("out_table.tsv")
  
#RBCK1 comes up twice for eosinophil # and %. 
ca_results_var_window %>% 
  group_by(conditional_locus, chr, phenotype, pheno_class) %>% 
  summarise(var_at_locus = sum(vars_in_window), genes_in_locus = n()) %>% 
  filter(var_at_locus == 0)  %>% 
  ungroup() %>% 
  count(conditional_locus)



### Variants within a Blood Cell Class. 

s <- ca_results_var_window_byClass %>%
  group_by(locus_name, chr, phenotype, pheno_class)  %>% 
  summarise(var_at_locus = sum(vars_in_window),
            var_at_locus_byClass = sum(vars_in_window_byClass), 
            genes_in_locus = n()) %>% 
  filter(var_at_locus_byClass == 0) %>% 
  select(-var_at_locus) %>% 
  ungroup()


sub <- ca_results_var_window_byClass %>%
  group_by(locus_name, chr, phenotype, pheno_class)  %>% 
  summarise(var_at_locus = sum(vars_in_window),
            var_at_locus_byClass = sum(vars_in_window_byClass), 
            genes_in_locus = n()) %>% 
  filter(var_at_locus_byClass > 0)


# 70 loci with no known variants in window of any genes at locus. 
ca_results_var_window_by_loci <- ca_results_var_window_byClass %>%
  group_by(conditional_locus, chr, phenotype, pheno_class) %>% 
  summarise(var_at_locus = sum(vars_in_window),
            var_at_locus_byClass = sum(vars_in_window_byClass), 
            genes_in_locus = n()) %>% 
  filter(var_at_locus_byClass == 0) %>% 
  select(-var_at_locus) %>% 
  ungroup()

# 57 unique genes that are sentinel's at these loci. 
ca_results_var_window_by_loci %>% 
  distinct(conditional_locus, locus_name, phenotype)

#13 loci with no known variants and multiple genes. 
cond_list <- ca_results_var_window_by_loci %>% 
  arrange(desc(genes_in_locus))  %>% 
  select(conditional_locus, phenotype) 

#Merge to go from loci level to gene level
ca_results_var_window_byClass_sub <- inner_join(ca_results_var_window_byClass, 
           cond_list)

ca_results_var_window_byClass_sub %>% 
  write_tsv("../../../data/functional_annotation_TWAS_results/novel_genes_92_with_no_GWAS_variants_by_class.tsv")

inner_join(ca_results_var_window_byClass, sub) %>% 
  filter(conditional_locus == "EPO",
         phenotype == "hemoglobin_concentration") %>% 
  rowwise() %>% 
  mutate(vars_in_window_data = list(pull_novel(.chr = chr,
                                          .start_pos = start_pos,
                                          .end_pos = end_pos))) %>% 
  select(conditional_locus, phenotype, locus_name, vars_in_window_data) %>% 
  unnest(cols = vars_in_window_data) %>% 
  mutate(var_reported_pheno_class = case_when(
    str_detect(var_reported_pheno,"baso|eos|lymph|neutr|mono|white") ~ "WBC", 
    str_detect(var_reported_pheno, "plate") ~ "PLT", 
    TRUE ~ "RBC"
  )) %>% 
  write_tsv("../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_novel_lociByClass_plotData_varsEPO.tsv")


write_tsv(x = ca_results_var_window_byClass_sub, 
          path = "../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_novel_lociByClass_plotData.tsv")

