library(tidyverse)


novel_genes <- read_tsv("../../../data/functional_annotation_TWAS_results/no_known_gwas_variants_genes.tsv")


marg_betas <- read_tsv("../../../data/analysis_results/regenie/association_results/marg_results_w_beta.tsv") %>% 
  select(gene_name, phenotype, marg_beta, 
         marg_se)

mvp <- read_tsv("../../../data/MVP_replication_results/mvp_rep_full.tsv") %>% 
  select(gene_name, phenotype, contains("mvp"))


get_pheno_betas <- function(gene, pheno){
  read_delim(sprintf("../../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_conditionalAnalysis_%s_results", 
                     pheno), 
             delim = " ") %>% 
    filter(ID == gene) %>% 
    select(conditional_beta = BETA, conditional_se = SE)
}

conditional_betas <- novel_genes %>% 
  distinct(gene_name, phenotype) %>% 
  rowwise() %>% 
  mutate(d = list(get_pheno_betas(gene_name, phenotype))) %>%
  unnest(d)

tss <- read_tsv("../../../../reference_files/gencode_v28_b38_tss", 
                col_names = c("chr", "tss", "strand", "gene_name")) %>% 
  select(gene_name, tss)


gene_table <- novel_genes %>% 
  left_join(marg_betas) %>% 
  left_join(conditional_betas) %>% 
  left_join(mvp) %>% 
  left_join(tss) %>% 
  mutate(log10_mvp_p = -log10(mvp_p), 
         in_mvp = if_else(is.na(mvp_beta), 0, 1),
         `TWAS Beta` = str_c(round(marg_beta, 3), " (", round(marg_se, 3), ")"), 
         `CA Beta` = str_c(round(conditional_beta, 3), " (", round(conditional_beta, 3), ")"), 
         `MVP Beta` = str_c(round(mvp_beta, 3), " (", round(mvp_se, 3), ")"), 
         phenotype = phenotype %>% 
           str_replace_all("_", " ") %>% 
           str_to_title()) %>% 
  arrange(desc(in_mvp), phenotype, chr, start_pos, end_pos) %>% 
  select(-c(in_mvp, marg_beta, marg_se, conditional_beta, conditional_se, mvp_beta, mvp_se)) %>% 
  select(phenotype, gene_name, chr, tss, model_r2, `TWAS Beta`, log10_regenie_p, `CA Beta`, log10_conditional_p, 
         `MVP Beta`, log10_mvp_p)%>% 
  rename_with(.fn = ~str_replace(.x, "_", " ")) %>% 
  rename_with(str_to_title)

write_tsv(gene_table, "../../../data/tables/novel_gene_table.tsv")
