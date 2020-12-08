#!/bin/bash

for chr in 1 6 19; do
	sbatch -J gene_var_cor_chr${chr} -t 12:00:00 --mem 250MB --mail-type END --mail-user bryce.rowland@unc.edu --out gene_var_cor_chr${chr}.out --wrap "../compute_var_gene_cor $chr"
done

# for chr in {2..5} {7..18} {20..22}; do
#         sbatch -J gene_var_cor_chr${chr} -t 4:00:00 --mem 250MB --mail-type END --mail-user bryce.rowland@unc.edu --out gene_var_cor_chr${chr}.out --wrap "../compute_var_gene_cor $chr"
# done
