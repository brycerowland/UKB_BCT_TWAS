## Training DGN Gene Expression Prediction Models Using an Elastic Net Pipeline.

The scripts in this directory are for pre-processing the imputed genotypes in DGN and UKB, training gene expression prediction models, and then using them to impute gene expression in UK Biobank. Scripts in this directory will be referred in this README by their numerical prefix. For example `01_generate_subset_maf_files_b38_chunks.sh` will be referred to via shorthand as `01`.

Scripts `00a`, `00b`, `01`, and `02` each perform variant matching and quality control for SNPs in DGN and UKB. To maximize variant overlap between DGN and UKB, we only train DGN models on variants that are QC'd in both DGN and UKB. 

Scripts `03` and `04` are helper scripts to get reference files into the correct format to run the gene expression training pipeline in script `05a`. Script `05a` executes `R/EN_0.8.setseed.R`, which is the script containing the details of our gene expression model training procedure. In summary, we train models to predict gene expression using SNPs within +/- 1Mb of a gene using an Elastic Net model with alpha = 0.5. 

Script `06` is used to compute CV R2 for the gene expression prediction models (results not included in our work, but could be useful for future researchers). Scripts `07`, `08`, and `09` each create files and then use files to impute gene expression in UKB using the DGN trained models. 
