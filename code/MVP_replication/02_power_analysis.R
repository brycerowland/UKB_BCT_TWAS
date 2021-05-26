library(pwr)
library(broom)

twas_power <- function(beta_hat, 
                       alpha){
  if(sign(beta_hat)==1){
    alt <- "greater"
  } else {
    alt <- "less"
  }
  pwr.t.test(n = 141286, 
             d = beta_hat, 
             sig.level =  alpha,
             alternative = alt, 
             type = "one.sample") %>% 
    tidy() %>% 
    pull(power)
}

int_data %>% arrange(gene_name) %>% 
  View()

power_df <- left_join(int_data, mvp_full, 
          by = c("gene_name", "phenotype")) %>% 
  mutate(mvp_sign_match = if_else(sign(mvp_beta) == sign(marg_beta), 
                                  1, 0)) 

marg_results %>% 
  filter(log10_regenie_p > 6.763558) %>% 
  count(abs(marg_beta) > 0.0137232)


marg_power_df <- left_join(new_marg, marg_results %>% 
            select(-log10_regenie_p), 
                           by=c("gene_name", "phenotype")) %>% 
  left_join(mvp_full) %>% 
  filter(log10_regenie_p > 6.763558) %>% 
  filter(pheno_not_in_mvp == 0, 
         is.na(twas_model_not_in_mvp)) %>% 
  mutate(mvp_sign_match = if_else(sign(mvp_beta) == sign(marg_beta), 
                                  1, 0)) %>% 
  mutate(r2_cat = cut(model_r2, breaks = c(0.05,0.1,0.2,0.4,0.6, 1)))



marg_twas_power <- marg_power_df %>% 
  rowwise() %>% 
  mutate(marg_power = twas_power(marg_beta, 
                                 alpha = 0.05/5993)) %>%
  ungroup() %>% 
  mutate(marg_rep_new_bon = if_else(mvp_p < 0.05/5993 & mvp_sign_match == 1, 
                                1, 0), 
         marg_rep_new_nom = if_else(mvp_p < 0.05 & mvp_sign_match == 1, 
                                                                  1, 0))

marg_twas_power %>% 
  summarise(rep_rate_bon = mean(marg_rep_new_bon), 
            rep_rate_nom = mean(marg_rep_new_nom),
            n = n())
  

marg_twas_power %>% 
  group_by(r2_cat) %>% 
  summarise(rep_rate_bon = mean(marg_rep_new_bon), 
            rep_rate_nom = mean(marg_rep_new_nom),
            n = n())



marg_twas_power %>% 
  select(marg_beta, marg_se,
         mvp_beta, mvp_se, mvp_p, 
         marg_rep_new)

pwr.t.test(n = 141286, sig.level = 0.05/5993, 
            type = "one.sample", 
           alternative = "greater", 
           power = .8)

marg_twas_power %>% 
  ggplot(aes(x = marg_beta)) + 
  geom_histogram() + 
  geom_vline(xintercept = 0.0137232, 
             color = "red", 
             linetype = "dashed") + 
  geom_vline(xintercept = -1 * 0.0137232, 
             color = "red", 
             linetype = "dashed")


marg_twas_power %>% 
  mutate(marg_rep_new = if_else(mvp_p < 0.05/5993 & mvp_sign_match == 1, 
                                1, 0)) %>% 
  ungroup() %>% 
  group_by(phenotype) %>% 
  count(marg_rep_new)  %>% 
  mutate(pct = n/sum(n)) %>% 
  filter(marg_rep_new == 1) %>% 
  arrange(pct)


marg_twas_power 
%>% 
  ggplot(aes(x = marg_beta, 
             y = mvp_beta,
             color = as.factor(marg_rep_new_nom))) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~r2_cat)

marg_twas_power %>% 
  select(marg_power) %>% 
  mutate(lt_0.2 = if_else(marg_power < 0.2, 1, 0), 
         lt_0.5 = if_else(marg_power < 0.5, 1, 0), 
         gt_0.5 = if_else(marg_power > 0.5, 1, 0), 
         gt_0.8 = if_else(marg_power > 0.8, 1, 0)) %>% 
  summarize(across(.cols = -marg_power, 
                  .fns = sum)) %>% 
  pivot_longer(cols = everything(),
               names_to = "power_cat", 
               values_to = "n") %>% 
  mutate(pct = n/5993)

marg_twas_power %>% 
  ggplot(aes(x = marg_beta)) + 
  geom_histogram() 


power_df %>% 
  filter(novel_gene == 1) %>% 
  select(gene_name, phenotype, marg_beta, marg_se, mvp_beta, 
         mvp_nominal_rep, mvp_sign_match)  %>% 
  rowwise() %>% 
  mutate(power = twas_power(marg_beta))

power_df%>% 
  rowwise() %>% 
  mutate(power = twas_power(marg_beta, 
                            alpha = 0.5/308), 
         mvp_rep_new = if_else(mvp_nominal_rep == 1 && mvp_sign_match == 1, 
                               "rep", "no rep")) %>% ungroup() %>% 
  filter(gene_name == "SNHG5") %>% 
  select(marg_beta, gene_name, phenotype, power)

power_df %>% 
  count(pheno_not_in_mvp, twas_model_not_in_mvp)


cond_power <- power_df %>% 
  filter(pheno_not_in_mvp == 0, 
         is.na(twas_model_not_in_mvp)) %>% 
  rowwise() %>% 
  mutate(power = twas_power(marg_beta, 
                            0.05/301), 
         mvp_rep_new_nom = if_else(mvp_nominal_rep == 1 && mvp_sign_match == 1, 
                               "rep", "no rep"), 
         mvp_rep_new_bon = if_else(mvp_rep_bonCA == 1 && mvp_sign_match == 1, 
                                   "rep", "no rep")) 
cond_power %>%
  ungroup() %>% 
  select(power)%>% 
  mutate(lt_0.2 = if_else(power < 0.2, 1, 0), 
         lt_0.5 = if_else(power < 0.5, 1, 0), 
         gt_0.5 = if_else(power > 0.5, 1, 0), 
         gt_0.8 = if_else(power > 0.8, 1, 0)) %>% 
  summarize(across(.cols = -power, 
                   .fns = sum)) %>% 
  pivot_longer(cols = everything(),
               names_to = "power_cat", 
               values_to = "n") %>% 
  mutate(pct = n/301)

marg_twas_power %>% 
  filter(gene_name == "SNHG5") %>% 
  select(gene_name, phenotype, marg_beta, marg_power)


cond_power %>% 
  filter(novel_gene == 1) %>% select(chr:marg_beta, power, 
                                     mvp_rep_new_nom, mvp_rep_new_bon, 
                                     mvp_p)



int_data %>% 
  filter(gene_name == "IRAK1BP1") %>% 
  View()

View(marg_results %>% 
       arrange(gene_name))

power_df$mvp_p

qqunif(power_df$mvp_p)

