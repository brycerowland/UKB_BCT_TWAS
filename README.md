# Analysis code for UKB TWAS of 29 Blood Cell Traits

Here, we provide the analysis scripts used in our recent TWAS of blood cell traits in UK Biobank Europeans. Currently, our work is available as a preprint here: https://www.biorxiv.org/content/10.1101/2021.08.03.453690v2.abstract. We thank you for your interest in this work!

## Methods Overview
Our TWAS consists of several key analytical steps. First, we trained gene expression prediction models using an Elastic Net pipeline using RNA-seq data from 922 Europeans in the DGN cohort (`code/genotype_preprocessing`). Second, we tested blood cell trait phenotypes (`code/phenotype_preprocessing`) with imputed gene expression in UKB Europeans using REGENIE (`code/REGENIE_marginal_association`). Third, we performed several rounds of secondary analysis including conditional analysis to discover novel TWAS loci (`code/conditional_analysis/REGENIE`), TWAS finemapping (`code/finemapping`), and pathway analysis (`code/pathway_analysis`). Finally, we systematically assigned GWAS loci from a large GWAS study of blood cell traits to our TWAS genes (`code/twas_gene_assignment_for_gwas_variants`).

Data cleaning scripts from our MVP replication analysis are contained in `code/MVP_replication`, and code to generate supplementary tables are in `code/tables`.

## Analysis Code

