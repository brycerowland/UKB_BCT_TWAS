library(tidyverse)
library(ggvenn)
library(ggsci)
library(cowplot)

twas_scores <- read_tsv("twas_pct_table.tsv")
coloc_scores <- read_tsv("coloc_pct_table.tsv")
blueprint_scores <- read_tsv("pct_blueprint_table.tsv")
venn_data <- readRDS("../../../data/finemapping/venn_data.rds")
names(venn_data) <- c("Coloc", "TWAS")


coloc_v2g <- read_tsv("../gwas_variants_coloc_gene_assignments.tsv")
coloc_phenos <- coloc_v2g %>% 
  pull(UKB_TWAS_trait) %>% 
  unique()

twas_v2g <- read_tsv("../gwas_variants_twas_gene_assignments.tsv") %>% 
  filter(twas_trait %in% coloc_phenos)


coloc_hist_df <- coloc_v2g %>% 
  group_by(UKB_TWAS_trait, condsig_var) %>% 
  count() %>% 
  rename(twas_trait = UKB_TWAS_trait, 
         unique_variant_id_b37 = condsig_var, 
         n_genes = n) %>% 
  mutate(method = "coloc")

twas_hist_df <- twas_v2g %>% 
  group_by(twas_trait, unique_variant_id_b37) %>% 
  count() %>% 
  rename(n_genes = n) %>% 
  mutate(method = "TWAS")

bind_rows(coloc_hist_df, 
          twas_hist_df) %>% 
  group_by(method) %>% 
  summarise(mean(n_genes), 
            sd(n_genes), 
            median(n_genes)) 

a <- bind_rows(coloc_hist_df, 
          twas_hist_df) %>% 
  mutate(n_genes = as.factor(n_genes)) %>% 
  ggplot(aes(x = n_genes)) + 
  geom_bar(aes(x = method, 
               fill = n_genes, 
               group = n_genes), 
           position = "fill") + 
  theme_bw() + 
  labs(x = "Method", 
       y = "Proportion of genes assigned
       per association") + 
  scale_fill_discrete(name = "Genes assigned\nper association")


pd <- bind_rows(twas_scores,
          coloc_scores, 
          blueprint_scores) %>% 
  mutate(ref = case_when(
    ref == "csq" ~ "Soranzo\nCSQ", 
    ref == "gene_exp" ~ "Soranzo\nExp.", 
    ref == "ot_any" ~ "OT\nAny", 
    ref == "ot_max" ~ "OT\nMax", 
    ref == "blueprint_exact_match" ~ "BLUEPRINT\nExact", 
    ref == "blueprint_some_match" ~ "BLUEPRINT\nAny"
  )) %>% 
  filter(!str_detect(ref, "Soranzo")) %>% 
  mutate(method = if_else(method == "twas", "TWAS", "coloc"), 
         method = fct_relevel(method, rev))




b <- pd %>%  
  ggplot(aes(x = ref,
             y = n, 
             fill = method)) + 
  geom_bar(stat = "identity", 
           position = "dodge")  + 
  theme_bw() + 
  labs(x = "Reference\nDataset", 
       y = "Count") + 
  scale_fill_npg(name = "Method")


c <-pd %>%  
  ggplot(aes(x = ref, 
             group = method, 
             fill = method)) + 
  geom_bar(aes(y = pct), 
           stat = "identity", 
           position = "dodge") + 
  theme_bw() + 
  labs(x = "Reference\nDataset", 
       y = "Percent") + 
  scale_fill_npg(name = "Method")



venn <- ggvenn(venn_data, 
               c("Coloc", "TWAS"),
               fill_color = rev(pal_npg("nrc")(2)), 
               set_name_size = 4,
               text_size = 2) + theme_bw() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

prow <- plot_grid(venn,
          b+ theme(legend.position="none", text = element_text(size = 8)),
          c+ theme(legend.position="none", text = element_text(size = 8)), 
          nrow = 1, 
          align = 'h', 
          labels = c("A", "B", "C"))



legend <- get_legend(b)
fig5 <- plot_grid(prow, legend, rel_widths = c(3, .4)) 
fig5
save_plot("figure5_v2g.jpg", 
          fig5, 
          base_width = 174, 
          base_height = 50,
          units = "mm", 
          dpi = 300)


