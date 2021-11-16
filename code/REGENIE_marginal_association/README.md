## REGENIE Association Testing. 

We highly recommend that interested users check out the excellent REGENIE documentation on their website before evaluating our association testing pipeline: https://rgcgithub.github.io/regenie/. 

### Code

Scripts `01`-`04` are concerned with preprocessing the necessary files to run REGENIE, including the imputed gene expression in UKB and hard call genotypes to control for population stratification. 

The association testing in REGENIE consists of two steps as outlined in their documentations, and are run in the `regenie_step1.sh` and `regenie_step2.sh` files. 

Script `05` assembles our TWAS results. Script `06` assigns each TWAS gene to a TWAS loci. In breif, each gene is assigned to a TWAS sentinel gene within +/- 1Mb, until all significant genes are in a TWAS locus. 
