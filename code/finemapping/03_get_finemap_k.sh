#!/bin/bash

for pheno in white_blood_cell_count red_blood_cell_count hemoglobin_concentration hematocrit_percentage mean_cell_volume mean_cell_hemoglobin mean_cell_hemoglobin_concentration red_blood_cell_width platelet_count mean_platelet_volume plateletcrit platelet_distribution_width lymphocyte_count monocyte_count neutrophil_count eosinophil_count basophil_count lymphocyte_percentage monocyte_percentage neutrophil_percentage eosinophil_percentage basophil_percentage reticulocyte_percentage reticulocyte_count mean_reticulocyte_volume mean_sphered_cell_volume immature_reticulocyte_fraction hls_reticulocyte_percentage hls_reticulocyte_count; do

echo $pheno

locus_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results
finemap_dir=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/FINEMAP/output_files/${pheno}

rm -f k_vector_${pheno}
rm -f prob_vector_${pheno}

while read locus_name n_genes; do
	
	if [ $n_genes -eq 1 ]; then
		echo NA >> k_vector_${pheno}
        	echo NA >> prob_vector_${pheno}
		continue	

	fi
	n_cred_files=$( ls ${finemap_dir}/*${locus_name}.cred* | wc -l )	
	
	#initialize best_k
	best_k=-1
	best_prob=-1



	for file in ${finemap_dir}/*${locus_name}.cred*; do
		
		if [ $n_cred_files -eq 1 ]; then
			best_k=$( basename $file | sed 's/.*.cred//g' )
			best_prob=$( head -1 $file | cut -f 2 -d "=" | sed 's/^ *//g' )
		else
			k=$( basename $file | sed 's/.*.cred//g' )
			prob=$( head -1 $file | cut -f 2 -d "=" | sed 's/^ *//g' )
			if (( $(echo "$prob > $best_prob" |bc -l) )); then
				best_k=$k
				best_prob=$prob
			fi
	
		fi	


	done
	echo $best_k >> k_vector_${pheno}
	echo $best_prob >> prob_vector_${pheno}

done < ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_locusInfo 

paste ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_locusInfo k_vector_${pheno} prob_vector_${pheno} > ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_locusInfo_FINEMAP

done

