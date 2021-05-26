library(tidyverse)


read_tsv("../twas_gene_assignment_for_gwas_variants/gwas_variants_twas_gene_assignments.tsv") %>% 
  select(unique_variant_id_b37,
         phenotype = twas_trait, 
         pheno_class = associated_blood_index_class, 
         gene_name = twas_gene, 
         r2) %>% 
  write_tsv("../../data/tables/variant_to_gene_supp_table.tsv")
