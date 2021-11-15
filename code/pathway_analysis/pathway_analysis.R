library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

finemap <- read_tsv("../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_FINEMAP_suppTable.tsv")


finemap %>% 
  filter(!is.na(prob)) %>% 
  filter(prob > 0.5)


marg_twas <- read_tsv("../../data/analysis_results/regenie/association_results/UKB_BCT_TWAS_results_FINEMAP_ref_file_TWASLoci_noHLA.tsv")

marg_w_beta <- read_tsv("../../data/analysis_results/regenie/association_results/marg_results_w_beta.tsv")

keytypes(org.Hs.eg.db)


twas_go <- function(twas_pheno){
  twas_sig <- marg_twas %>% 
    filter(sig_high_PPCS == 1, 
           phenotype == twas_pheno)

  
  sig_w_beta <- inner_join(marg_w_beta, twas_sig)
  
  universe_genes <- marg_w_beta %>% distinct(gene_name) %>% pull(gene_name)
  
  new_sig <- bitr(sig_w_beta$gene_name, 
                  fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  new_universe <- bitr(universe_genes, 
                       fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  
  egobp <- clusterProfiler::enrichGO(
    gene     = new_sig[[2]],
    universe = new_universe[[2]],
    OrgDb    = org.Hs.eg.db,
    ont      = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05,
    minGSSize = 10,
    readable = TRUE)
  
  if(egobp@result %>% filter(p.adjust < 0.05) %>% nrow() > 0){
    print("enriched GO terms")
    return(egobp)
  } else{
    print("no enriched GO terms")
    return(NULL)
  }
}


#Run pathway analysis for each phenotype. 
pathway_results <- marg_twas %>%
  distinct(phenotype) %>% 
  rowwise() %>% 
  mutate(pathway_results = list(twas_go(phenotype)))

  

sig_terms <- pathway_results %>% 
  filter(!is.null(pathway_results)) %>% 
  mutate(sig_terms = list(pathway_results@result %>% 
           filter(p.adjust < 0.05) %>% 
           as_tibble()), 
         n_go_terms = nrow(sig_terms)) %>% 
  arrange(desc(n_go_terms)) 

sig_terms %>% 
  unnest(sig_terms) %>% 
  select(Phenotype = phenotype, ID, Description, `Gene Ratio`=GeneRatio, `Background Ratio` = BgRatio, 
         `FDR p-value` = p.adjust, `GO Term Genes` = geneID) %>% 
  write_tsv("../../data/tables/pathway_analysis_sup_table.tsv")






get_goplots <- function(out_path, 
                        results, pheno, n_go_terms){
  p <- results %>% 
    filter(phenotype == pheno) %>% 
    pull(pathway_results) %>% 
    .[[1]] %>% 
    goplot(showCategory = n_go_terms, 
           font.size = 0.5) + theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_blank(), 
          axis.ticks = element_blank())

  save_plot(plot = p,
         filename = sprintf("%s/%s_pathway_results_nGO_%s.jpg", 
                            out_path, pheno, n_go_terms), 
         base_width = 6, 
         base_height = 6 * 9 / 16)
}


twas_phenos <- sig_terms %>% filter(n_go_terms > 1) %>% distinct(phenotype) %>% pull(phenotype)

map(twas_phenos, ~get_goplots("../../data/pathway_analysis/plots", 
                              pathway_results, .x, n_go_terms = 5))

get_goplots("../../data/pathway_analysis/plots/", 
            pathway_results, "lymphocyte_count", n_go_terms = 5)


