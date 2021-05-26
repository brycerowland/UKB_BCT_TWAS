#!/usr/bin/env Rscript

library(tidyverse)
library(janitor)

group_chr <- function(data, chr_num, window){
  chr <- data %>% 
    filter(chr == chr_num) %>% 
    arrange(start_pos) %>% 
    mutate(locus_name = "")
  
  sig_groups <- chr %>% filter(log10_regenie_p > 6.763558) %>% pull(locus_name)
  
  i <- 0
  #Loop until all TWAS significant genes are assigned a locus. 
  while(any(sig_groups == "")){
    #loop at top TWAS p-value for gene not in a group 
    top_gene <- chr %>% filter(locus_name=="") %>% arrange(desc(log10_regenie_p)) %>%
      slice(1) %>%
      pull(gene_name)
    
    #get positions
    s <- chr %>% filter(gene_name == top_gene) %>% pull(start_pos)
    e <- chr %>% filter(gene_name == top_gene) %>% pull(end_pos)
    
    #Get genes in group for genes not yet in a group
    gene_group <- chr %>% 
      filter(start_pos > s - window,
             end_pos < e + window, 
             locus_name == "") %>% 
      pull(gene_name)
    
    chr <- chr %>% 
      mutate(locus_name = case_when(
        gene_name %in% gene_group ~ top_gene, 
        TRUE ~ locus_name
      ))
    sig_groups <- chr %>% filter(log10_regenie_p > 6.763558) %>% pull(locus_name)
    i <- i + 1
  }
  chr <- chr %>% 
    mutate(locus_name = if_else(locus_name == "", NA_character_, locus_name))
  return(chr)
}

group_genome <- function(data, window){
  chr_list <- data %>% 
    pull(chr) %>% 
    unique()
  
  l <- vector("list", length = length(chr_list))
  
  for(i in chr_list){
    l[[i]] <- group_chr(data, i, window)
  }
  
  
  return(bind_rows(l))
}

count_groups <- function(data){
  data %>% 
    pull(locus_name) %>% 
    unique() %>% 
    length()
}

group_phenos <- function(pheno){
  case_when(
    pheno %in% c("basophil_count", "eosinophil_count",
                 "lymphocyte_count", "monocyte_count", "neutrophil_count",
                 "white_blood_cell_count") ~ "WBC",
    pheno %in% c("platelet_count", "mean_platelet_volume") ~ "PLT",
    TRUE ~ "RBC"
  )
}

get_sig_genes <- function(data){
  data %>% 
    pull(gene_name)%>% 
    unique()
}

get_lead_genes <- function(data){
  data %>% 
    pull(locus_name)%>% 
    unique()
}


files <- list.files(path = "../../data/analysis_results/regenie/association_results", pattern = "_ref_file", 
           recursive = T,
           full.names = T)
files <- files[!grepl("bonferroniSigGenes", files)]


for(file_path in files){
  print(file_path)
  results <- read_tsv(file = file_path) %>% 
    clean_names()
  write_tsv(path = str_c(file_path, "_TWASLoci"), 
            x = group_genome(results, 1e6))
}






