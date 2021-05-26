#!/usr/bin/env Rscript

library(tidyverse)

pheno <- read_tsv("/pine/scr/b/r/bryce38/twas/UKB_temp/data/phenotype_data/UKB_BloodCellTraits_EUR_TWAS_w_varNames.tsv")
plot_dir <- "/proj/yunligrp/users/bryce/twas/UKB/data/phenotype_data/pheno_plots"


hist_fxn <- function(data, var){
  file_name <- str_c(plot_dir,"/", var, "_histogram.png")
  ggplot(data, aes(x = .data[[var]])) + geom_histogram() +
    labs(x = var)
  ggsave(filename = file_name)
}

map(colnames(pheno)[c(45:54,56:74)], ~hist_fxn(pheno, .x))
