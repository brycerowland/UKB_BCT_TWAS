library(tidyverse)
library(janitor)
library(cowplot)
library(ggrepel)

files <- list.files("../../../data/conditional_analysis/REGENIE/analysis_results/",
                    pattern = "_refFile$", 
                    full.names = T)


l <- vector("list", length = length(files))

for(i in 1:length(files)){

  l[[i]] <- read_delim(files[[i]], delim= ' ') %>% 
    clean_names()

}

#$1==6 && $2 > 24477797 && $3 < 37448354

conditional_results <- bind_rows(l) %>% 
  filter(!(chr==6 & start_pos > 24477797 & end_pos < 37448354)) %>% 
  mutate(pheno_class = case_when(
    phenotype %in% c("basophil_count", "basophil_percentage",
                     "eosinophil_count", "eosinophil_percentage", 
                     "lymphocyte_count", "lymphocyte_percentage",
                     "monocyte_count", "monocyte_percentage", 
                     "neutrophil_count", "neutrophil_percentage", 
                     "white_blood_cell_count") ~ "WBC", 
    phenotype %in% c("platelet_count", "mean_platelet_volume", 
                     "platelet_distribution_width", "plateletcrit") ~ "PLT", 
    TRUE ~ "RBC"
  ))

lag_vals <- conditional_results %>% 
  select(chr, start_pos) %>% 
  group_by(chr) %>% 
  filter(start_pos == max(start_pos)) %>% 
  distinct() %>% 
  arrange(chr) %>% 
  ungroup() %>% 
  mutate(lag_max = cumsum(lag(start_pos, default = 0)))  %>% 
  select(chr, lag_max)

cond_results_cum <- left_join(conditional_results, lag_vals, by = "chr") %>% 
  mutate(cum_pos = start_pos + lag_max) 

axis_ticks <- cond_results_cum %>% 
  group_by(chr) %>% 
  summarize(pos = (max(cum_pos) + min(cum_pos))/2) %>% 
  pull(pos)

chr <- 1:22
chr[chr %% 2 == 0] <- ""
axis_labels <- chr %>% 
  as.character()

plot_data <- cond_results_cum %>% 
  mutate(is_even = if_else((chr %% 2) == 0, 
                           1, 0)) %>% 
  group_by(gene_name, pheno_class) %>% 
  mutate(min_p = max(log10_conditional_p)) %>% 
  filter(log10_conditional_p == min_p) %>% 
  ungroup() 

# write_tsv(path = "../../../data/conditional_analysis/plots/conditional_analysis_plot_data.tsv", 
          # plot_data)

font_size_theme <- theme(text = element_text(size = 6))

cond_plot <- plot_data %>% 
  ggplot(aes(x = cum_pos, y = log10_conditional_p, 
             color = is_even)) + 
  geom_text_repel(aes(label = if_else(log10_conditional_p > 12, gene_name, "")), 
                  size = 2) + 
  geom_point(alpha = 0.75, 
             size = 1) + 
  guides(color = F) + 
  facet_wrap(~pheno_class) + 
  theme_bw() + 
  labs(y = "-log10 Conditional p-value", 
       x = "Chromosome", 
       title = "") + 
  scale_x_continuous(breaks = axis_ticks, 
                     labels = axis_labels)+ 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 5.434152, 
             linetype = "dashed", 
             color = "red") + 
  font_size_theme
cond_plot



ggsave(filename = "../../../data/conditional_analysis/plots/fig3_man_plot_CA.tiff",
          plot = cond_plot, 
          width = 174, 
          height = 174 * (3/5),
       dpi = 1000,
       units = "mm")


pacsin2 <- cond_results_cum %>% 
  filter(pheno_class == "PLT", 
         locus_name == "PACSIN2", 
         phenotype == "platelet_distribution_width") %>% 
  select(gene_name, cum_pos, log10_regenie_p, 
         log10_conditional_p, phenotype, locus_name) %>% 
  pivot_longer(cols = starts_with("log"), 
              names_to ="test_type", 
              values_to = "p")



 marginal <- pacsin2 %>% 
  filter(test_type == "log10_regenie_p") %>% 
  mutate(label = if_else(p > 400, gene_name, "")) %>% 
  ggplot(aes(x = cum_pos, y = p)) + 
  geom_point() + 
  geom_text(aes(label = label), 
            nudge_x = 1e5) + theme_bw() + 
  labs(y = "marginal p", 
       x = "")

conditional <- pacsin2 %>% 
  filter(test_type == "log10_conditional_p")%>% 
  mutate(label = if_else(p > 30, gene_name, "")) %>% 
  ggplot(aes(x = cum_pos, y = p)) + 
  geom_point()  + 
  geom_text(aes(label = label), 
            nudge_x = 1e5) + theme_bw() + 
  scale_y_reverse() + labs(y = "conditional p",
                           x = "")


plot_grid(marginal, conditional, nrow = 2)




