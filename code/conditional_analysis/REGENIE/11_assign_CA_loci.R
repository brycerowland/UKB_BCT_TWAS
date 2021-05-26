#!/usr/bin/env Rscript

library(tidyverse)
library(janitor)

group_chr <- function(data, chr_num, window){
  chr <- data %>% 
    filter(chr == chr_num) %>% 
    arrange(start_pos) %>% 
    mutate(conditional_locus = "")
  
  sig_groups <- chr %>% filter(log10_conditional_p > 5.3714, 
                               log10_regenie_p > 6.763558) %>% pull(conditional_locus)
  
  i <- 0
  #Loop until all TWAS significant genes are assigned a locus. 
  while(any(sig_groups == "")){
    #loop at top TWAS p-value for gene not in a group 
    top_gene <- chr %>% filter(conditional_locus=="",
                               log10_regenie_p > 6.763558) %>% arrange(desc(log10_conditional_p)) %>%
      slice(1) %>%
      pull(gene_name)
    
    #get positions
    s <- chr %>% filter(gene_name == top_gene) %>% pull(start_pos)
    e <- chr %>% filter(gene_name == top_gene) %>% pull(end_pos)
    
    #Get genes in group for genes not yet in a group
    gene_group <- chr %>% 
      filter(((start_pos > (s - window) &
             end_pos < (e + window)) | 
               (start_pos < (s - window) &
                  end_pos > (s - window)) | 
               (start_pos < (e + window) & 
                  end_pos > (e + window))), 
             conditional_locus == "") %>% 
      pull(gene_name)
    
    chr <- chr %>% 
      mutate(conditional_locus = case_when(
        gene_name %in% gene_group ~ top_gene, 
        TRUE ~ conditional_locus
      ))
    sig_groups <- chr %>% filter(log10_conditional_p > 5.3714, 
                                 log10_regenie_p > 6.763558) %>% pull(conditional_locus)
    i <- i + 1
  }
  chr <- chr %>% 
    mutate(conditional_locus = if_else(conditional_locus == "", NA_character_, conditional_locus))
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
    pull(conditional_locus) %>% 
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
    pull(conditional_locus)%>% 
    unique()
}


files <- list.files(path = "../../../data/conditional_analysis/REGENIE/analysis_results/", pattern = "_refFile$", 
           recursive = T,
           full.names = T)



for(file_path in files){
  print(file_path)
  results <- read_delim(file = file_path, delim = " ") %>% 
    clean_names()
  write_tsv(file = str_c(file_path, "_CALoci"),
            x = group_genome(results, 1e6))

}




