library(tidyverse)
library(xlsx)



marg_beta <- read_tsv("../../data/analysis_results/regenie/association_results/marg_results_w_beta.tsv")
marg_results <- read_tsv("../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_FINEMAP_ref_file_TWASLoci_noHLA.tsv") %>% select(chr:gene_name, 
                                                                                                                                                        model_r2, phenotype,
                                                                                                                                                        locus_name:marginal_significant)

marg_finemap <- read_tsv("../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_FINEMAP_suppTable.tsv")  %>% 
  select(gene_name, phenotype, best_k)

marg_w_finemap <- left_join(marg_results,
          marg_beta, 
          by = c("gene_name", "phenotype")) %>% 
  left_join(marg_finemap, by = c("gene_name", "phenotype")) %>% 
  select(chr:phenotype, 
         marg_beta:log10_regenie_p, 
         marginal_significant,
         locus_name, best_k, prob) 

ca_results <- read_tsv('../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_CA_results.tsv') %>% 
  select(gene_name, phenotype, log10_conditional_p, conditional_locus)

pheno_key <- read_tsv('../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_CA_results.tsv') %>% 
  distinct(phenotype, pheno_class)


mvp_full <- read_tsv("../../data/MVP_replication_results/mvp_rep_full.tsv") %>% 
  select(-chr)


t <- left_join(marg_w_finemap, ca_results, 
          by = c("gene_name", "phenotype")) %>% 
  left_join(mvp_full, by = c("gene_name", "phenotype")) %>% 
  left_join(pheno_key) %>% 
  relocate(pheno_class, .before = phenotype) %>% 
  arrange(pheno_class, phenotype, chr, start_pos) %>% 
  filter(!is.na(locus_name)) 


t %>% 
  mutate(sent_gene = if_else(gene_name == locus_name & !is.na(prob), 1, 0)) %>% 
  group_by(locus_name, phenotype) %>% 
  summarise(s = sum(sent_gene)) %>% 
  filter(s == 0)

t %>% 
  distinct(locus_name, phenotype)


write_excel_csv(t, "../../data/tables/marginal_results_supp_table.csv")
