## Code

This subdirectory contains scripts that execute our TWAS conditional analysis. In breif, these scripts test the hypothesis that a TWAS signal is not driven by known GWAS variants within 1Mb of the TWAS significant gene. 

Scripts `01`-`06` perform preprocessing to construct the covariate file for each TWAS gene. For each TWAS significant gene-trait association, all sentinel variants from the Vuckovic et al GWAS for a blood cell trait within the same phenotype category are pulled as covariates. We then check their imputation quality and minor allele frequency in our analysis sample (`03`). We then construct an LD pruned set of variants (`04`), and extract variant dosage (`05`) in covariate format `06`. 

Conditional association testing is conducted in `07`-`09`, and results are processed in `10`-`12`. 
