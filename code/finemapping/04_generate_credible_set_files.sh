#!/bin/bash

for pheno in white_blood_cell_count red_blood_cell_count hemoglobin_concentration hematocrit_percentage mean_cell_volume mean_cell_hemoglobin mean_cell_hemoglobin_concentration red_blood_cell_width platelet_count mean_platelet_volume plateletcrit platelet_distribution_width lymphocyte_count monocyte_count neutrophil_count eosinophil_count basophil_count lymphocyte_percentage monocyte_percentage neutrophil_percentage eosinophil_percentage basophil_percentage reticulocyte_percentage reticulocyte_count mean_reticulocyte_volume mean_sphered_cell_volume immature_reticulocyte_fraction hls_reticulocyte_percentage hls_reticulocyte_count; do

finemap_dir=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/FINEMAP/output_files/${pheno}
locus_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/FINEMAP/best_k_credible_set_files/${pheno}

mkdir -p $out_dir

while read locus_name k; do

	#Process the best k credible set file.
	field_max=$(( 2 * $k ))
	for f in `seq 2 2 $field_max`; do 
		grep -v "^#" ${finemap_dir}/${pheno}_${locus_name}.cred${k} | tail -n +2 | cut -f $f,$(($f + 1)) -d ' '
	done | awk '$1!="NA"' > ${out_dir}/${pheno}_${locus_name}_${k}_causalGeneCredibleSet

	awk '$2 > .5' ${out_dir}/${pheno}_${locus_name}_${k}_causalGeneCredibleSet > ${out_dir}/${pheno}_${locus_name}_${k}_causalGeneCredibleSet_posteriorProb0.5

done < <( cat ${locus_dir}/UKB_BCT_TWAS_results_${pheno}_locusInfo_FINEMAP | awk '$2!=1' OFS="\t" | cut -f 1,3 )

done
