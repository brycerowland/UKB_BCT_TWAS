library(tidyverse)


fm_files <- list.files("~/bryce_group/twas/UKB/data/finemapping/FINEMAP/best_k_credible_set_files/", 
           pattern = "_causalGeneCredibleSet$", 
           recursive = T, 
           full.names = T)


l <- vector("list", 
            length = length(fm_files))


i <- 1

fm_files[i] %>% 
  basename() %>% 
  str_remove("_[0-9]_causalGeneCredibleSet") %>% 
  str_remove


read_delim(fm_files[i], delim = " ", 
           col_names = c("gene", "PPI"))

for(i in 1:length(fm_files)){
  locus_name <- fm_files[i] %>% 
    basename() %>% 
    str_remove("_[0-9]_causalGeneCredibleSet") %>% 
    str_remove(".*_")
  
  l[[i]] <- read_delim(fm_files[i], delim = " ", 
                       col_names = c("gene", "PPI"), 
                       col_types = cols(
                         gene = col_character(),
                         PPI = col_double()
                       )) %>% 
     mutate(locus_name = locus_name)
}


monica <- bind_rows(l)

monica
