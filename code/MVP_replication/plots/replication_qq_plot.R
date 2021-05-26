library(tidyverse)


rep <- read_tsv("../../../data/MVP_replication_results/ca_and_replication_results_allMargSigGenes.tsv")

mvp_qq_plot <- function(df, var, title){
  n <- nrow(df)
  expected <- tibble(
    expected = -log10(seq(from = 1, to  = 1/n, length.out = n))
  ) %>% 
    arrange(expected) %>% 
    mutate(rank = row_number())
  
  
  nonzero_min <- df %>% 
    filter({{var}} > 0) %>% 
    pull( {{ var }} ) %>% 
    min() 
  
  observed <- df %>% 
    filter(!is.na({{ var }})) %>% 
    mutate(observed = if_else({{var}} == 0, -log10(nonzero_min), -1 * log10({{var}}))) %>% 
    arrange(observed) %>% 
    mutate(rank = row_number())
  
  # browser()
  
  left_join(expected, observed, by = "rank") %>% 
    ggplot(aes(x = expected, y = observed)) + 
    geom_point() + 
    theme_bw() + 
    labs(x = "Expected -log10(p)", 
         y = "Observed -log10(p)", 
         title = title)
}



rep %>% 
  filter(!is.na(mvp_p)) %>% 
  mvp_qq_plot(mvp_p, 
              title  = "MVP Replication - TWAS marginally significant genes") 

rep %>% 
  filter(!is.na(mvp_p), 
         ca_sig == 1) %>% 
  mvp_qq_plot(mvp_p, 
              title = "MVP Replication - TWAS conditionally significant genes")


