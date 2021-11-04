library(tidyverse)
library(janitor)

rep_dir <- "../../data/MVP_replication_results/201226_fromGenisis_boltLMM_SNP_Results/"

rep_files <- list.files(path = rep_dir, 
           pattern = "*_res.txt", 
           full.names = T)
l <- vector("list", length = length(rep_files))

for(i in 1:length(rep_files)){
  mvp_pheno <- rep_files[[i]] %>% 
    basename() %>% 
    str_remove("boltLMMresultsSNPs_chrAll_") %>% 
    str_remove("_res.txt")
  
  l[[i]] <- read_tsv(rep_files[[i]], 
           col_types = cols(
             SNP = col_character(),
             CHR = col_double(),
             BP = col_double(),
             GENPOS = col_double(),
             ALLELE1 = col_character(),
             ALLELE0 = col_logical(),
             A1FREQ = col_double(),
             F_MISS = col_double(),
             BETA = col_double(),
             SE = col_double(),
             P_BOLT_LMM_INF = col_double(),
             P_BOLT_LMM = col_double()
           )) %>% 
    mutate(mvp_pheno = mvp_pheno)
}

df_rep_raw <- bind_rows(l) %>% 
  #Exclude HLA. Already was excluded. 
  filter(!(CHR == 6 & 
             GENPOS > 24477797 &
             GENPOS < 37448354))
df_rep_raw %>% 
  filter(SNP == "MFAP3L")


#Read in TWAS conditional analysis results 
ca_results <- read_tsv("../../data/conditional_analysis/REGENIE/analysis_results/UKB_BCT_TWAS_CA_results.tsv")%>% 
  mutate(ca_sig = if_else(log10_conditional_p > 5.3714, 1, 0))

df_rep <- df_rep_raw  %>% 
  mutate(phenotype = case_when(
    mvp_pheno == "Baso" ~ "basophil_count", 
    mvp_pheno == "Eos" ~ "eosinophil_count", 
    mvp_pheno == "HCT" ~ "hematocrit_percentage", 
    mvp_pheno == "Hemoglobin" ~ "hemoglobin_concentration", 
    mvp_pheno == "Lymph" ~ "lymphocyte_count", 
    mvp_pheno == "MCH" ~ "mean_cell_hemoglobin", 
    mvp_pheno == "MCHC" ~ "mean_cell_hemoglobin_concentration", 
    mvp_pheno == "MCV" ~ "mean_cell_volume", 
    mvp_pheno == "Mono" ~ "monocyte_count", 
    mvp_pheno == "MPV" ~ "mean_platelet_volume", 
    mvp_pheno == "Neut" ~ "neutrophil_count", 
    mvp_pheno == "Platelet" ~ "platelet_count", 
    mvp_pheno == "RBC" ~ "red_blood_cell_count", 
    mvp_pheno == "RDW" ~ "red_blood_cell_width", 
    mvp_pheno == "WBC" ~ "white_blood_cell_count"
  ), 
  pheno_class = case_when(
    str_detect(phenotype,"baso|eos|lymph|neutr|mono|white") ~ "WBC", 
    str_detect(phenotype, "plate") ~ "PLT", 
    TRUE ~ "RBC"
  )) %>% 
  select(gene_name = SNP, phenotype, pheno_class, 
         CHR, mvp_beta = BETA,
         mvp_se = SE, mvp_p = P_BOLT_LMM) %>% 
  clean_names()

write_tsv(df_rep, "../../data/MVP_replication_results/mvp_rep_full.tsv")

read_tsv("../../data/MVP_replication_results/mvp_rep_full.tsv")

mvp_phenos <- df_rep %>% distinct(phenotype) %>% pull()

# 5,993 TWAS significant genes for the 15 phenotypes for which we have matching
# data in MVP. 
# There are 349 genes which we were not able to replicate in MVP, even though we 
# had the phenotype. Follow up on these. 
rep_results <- left_join(ca_results, 
          df_rep, by = c("gene_name", 
                         "phenotype", "pheno_class",
                         "chr")) %>% 
  filter(phenotype %in% mvp_phenos) %>% 
  filter(!is.na(mvp_p)) %>% 
  mutate(mvp_rep = if_else(mvp_p < 0.05/557, 1, 0), 
         ca_sig = if_else(log10_conditional_p > 5.434152, 1, 0)) 




ca_results_w_rep_temp <- left_join(ca_results, 
                              df_rep, by = c("gene_name", 
                                             "phenotype", "chr")) %>% 
  mutate(mvp_rep_bonMarginal = if_else(mvp_p < 0.05/5993, 1, 0),
         mvp_rep_bonCA = if_else(mvp_p < 0.05/301, 1, 0), 
         mvp_nominal_rep = if_else(mvp_p < 0.05, 1, 0),
         pheno_not_in_mvp = if_else(!(phenotype %in% mvp_phenos), 
                                    1, 0)) 

missing_twas_model <- ca_results_w_rep_temp %>% 
  filter(pheno_not_in_mvp == 0, 
         is.na(mvp_p)) %>% 
  distinct(gene_name) %>% 
  mutate(twas_model_not_in_mvp = 1)

ca_results_w_rep <- left_join(ca_results_w_rep_temp, missing_twas_model)


# Out of 11,759 associations, 6,342 marginally significant associations have phenotypes in MVP
# 349 genes are missing predicition models from MVP.
# Thus, 5,993 marginally siginficant genes are eligable for replication. 
ca_results_w_rep%>% 
  filter(phenotype %in% mvp_phenos) %>% 
  filter(!is.na(mvp_p))

#Compute replication rates. 
ca_results_w_rep %>% 
  summarise(mvp_rep_total = sum(mvp_rep, na.rm = T), 
            mvp_nominal_rep_total = sum(mvp_nominal_rep, na.rm = T),
            mvp_rep_rate = mean(mvp_rep, na.rm = T), 
            mvp_nominal_rep_rate = mean(mvp_nominal_rep, na.rm = T))

write_tsv(ca_results_w_rep,"../../data/MVP_replication_results/ca_and_replication_results_allMargSigGenes.tsv")


#301 conditional sig genes available. Why not 557? 
## 335 have available phenotypes in MVP. 
## 34 are missing genes from MVP. 
#
## Of the 301 genes - 109 replicate in a strict sense (0.05/301) (36.2%).
# 203 in a nominal sense (67.4%)
ca_genes_in_mvp <- left_join(ca_results, 
          df_rep, by = c("gene_name", 
                         "phenotype", "chr"))  %>%
  filter(ca_sig == 1) %>% 
  filter(phenotype %in% mvp_phenos) %>% 
  filter(!is.na(mvp_p))

rep_results <- ca_genes_in_mvp %>% 
  mutate(mvp_rep = if_else(mvp_p < 0.05/nrow(.), 1, 0), 
         mvp_nominal_rep = if_else(mvp_p < 0.05, 1, 0)) 






#What about replications in any phenotypes in a pheno class? 
# There are 423 distinct genes*pheno_class
# 387 are available in MVP (1882 gene*pheno combinations).

ca_genes_by_class <- ca_results %>% 
  filter(ca_sig == 1) %>% 
  distinct(gene_name, pheno_class)


# 225/387 replicate at 0.005 level (multiple testing is bad here). 
# 115 replicate at 0.05/1882 level (1882 chances to find a "winner"). 
inner_join(df_rep, 
           ca_genes_by_class, 
           by = c("gene_name", "pheno_class")) %>% 
  group_by(gene_name, pheno_class) %>% 
  summarise(rep = sum(mvp_p < 0.05/1882), 
            rep_nominal = sum(mvp_p < 0.005)) %>% 
  filter(rep_nominal > 0) %>% 
  ungroup() %>% 
  summarise(sum(rep), 
            sum(rep_nominal))

################################# 
# What about our 9 novel signals??
################################# 

#Read in info for 9 "novel" genes
novel_genes <- read_csv("../conditional_analysis/REGENIE/copy_novel_genes.csv")

novel_genes %>% filter(!(phenotype %in% mvp_phenos))

#6 genes have a matching phenotype in MVP, although RBCK1 has both EOS# and EOS%.
novel_genes_in_mvp <- novel_genes %>% 
  distinct(chr, gene_name, phenotype) %>% 
  filter(phenotype %in% mvp_phenos)

#IRAK1BP1 is only gene with strict replication. 
inner_join(rep_results, 
           novel_genes_in_mvp, 
           by = c("chr", "gene_name", "phenotype")) %>% 
  select(phenotype, chr:log10_regenie_p, mvp_p, mvp_rep, 
         mvp_nominal_rep)

#No nominal replication for TMEM144
df_rep %>% filter(gene_name == "TMEM144")

#GTF2H2, which was conditionally significant for 
# platelet_distribution_width, replicates in MVP
# nominally for platelet_count.
df_rep %>% filter(gene_name == "GTF2H2")

#Weak statistical evidence for platelet count. 
df_rep %>% filter(gene_name == "MFAP3L")

#Weak statistical evidence for platelet_count. 
df_rep %>% filter(gene_name == "LPCAT4")

####################################################
# 82 with no nearby variants in related phenotypes #
####################################################

novel_genes_by_class <- read_tsv("../conditional_analysis/REGENIE/novel_genes_82_with_no_GWAS_variants_by_class.tsv")

#52 of these have phenotypes which strictly match to MVP. 
novel_genes_by_class_mvp <- novel_genes_by_class %>% 
  distinct(chr, gene_name, phenotype) %>% 
  filter(phenotype %in% mvp_phenos)

# 44/52 are in MVP. Why are 8 unavailable? 
# 16 out of the 44 have strict replication. 
# 32 out of 44 replicate at p = 0.05. 
inner_join(rep_results, 
           novel_genes_by_class_mvp, 
           by = c("chr", "gene_name", "phenotype")) %>% 
  filter(mvp_p < 0.05) %>% 
  select(phenotype, chr:log10_regenie_p, mvp_p, mvp_rep) %>% 
  print(n = 100)





