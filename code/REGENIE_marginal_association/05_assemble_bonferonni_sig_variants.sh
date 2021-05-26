#!/bin/bash

results_dir=/proj/yunligrp/users/bryce/twas/UKB/data/analysis_results/regenie/association_results
gene_ref_file=/proj/yunligrp/users/bryce/twas/UKB/data/prediction_files/gene_ref_files/gencode_build38_allChr_TWAS_genes.txt

for pheno in basophil_count basophil_percentage eosinophil_count eosinophil_percentage hematocrit_percentage hemoglobin_concentration hls_reticulocyte_count hls_reticulocyte_percentage immature_reticulocyte_fraction lymphocyte_count lymphocyte_percentage mean_cell_hemoglobin_concentration mean_cell_hemoglobin mean_cell_volume mean_platelet_volume mean_reticulocyte_volume mean_sphered_cell_volume monocyte_count monocyte_percentage neutrophil_count neutrophil_percentage platelet_count plateletcrit platelet_distribution_width red_blood_cell_count red_blood_cell_width reticulocyte_count reticulocyte_percentage white_blood_cell_count; do 
	echo $pheno

	#Create file of only Bonferonni sig genes. 
	cat <( head -1 ${results_dir}/UKB_BCT_TWAS_results_${pheno}.regenie | sed 's/ID/GENE/g' ) <( tail -n +2 ${results_dir}/UKB_BCT_TWAS_results_${pheno}.regenie | awk '$11 > 6.763558' ) > ${results_dir}/UKB_BCT_TWAS_results_${pheno}_bonferroniSigGenes

	#Add regenie p-value to gene ref file
	awk 'NR==FNR { a[$3]=$11;next } $4 in a { print $0, a[$4] }' OFS="\t"  ${results_dir}/UKB_BCT_TWAS_results_${pheno}_bonferroniSigGenes $gene_ref_file | sed '1 i\chr\tstart_pos\tend_pos\tgene_name\tDGN_chunk\tEN_fit\tmodel_r2\tcv_r2\tlog10_regenie_p'  > ${results_dir}/UKB_BCT_TWAS_results_${pheno}_bonferroniSigGenes_ref_file

	awk 'NR==FNR { a[$3]=$11;next } $4 in a { print $0, a[$4] }' OFS="\t"  ${results_dir}/UKB_BCT_TWAS_results_${pheno}.regenie $gene_ref_file | sed '1 i\chr\tstart_pos\tend_pos\tgene_name\tDGN_chunk\tEN_fit\tmodel_r2\tcv_r2\tlog10_regenie_p'  > ${results_dir}/UKB_BCT_TWAS_results_${pheno}_ref_file

	#Make phenotype category specific files. 
	if [ "$pheno" = "white_blood_cell_count" ] || [ "$pheno" = "lymphocyte_count" ] || [ "$pheno" = "monocyte_count" ] || [ "$pheno" = "neutrophil_count" ] || [ "$pheno" = "eosinophil_count" ] || [ "$pheno" = basophil_count ] || [ "$pheno" = "basophil_percentage" ] || [ "$pheno" = "eosinophil_percentage" ] || [ "$pheno" = "lymphocyte_percentage" ] || [ "$pheno" = "monocyte_percentage" ] || [ "$pheno" = "neutrophil_percentage" ];then
		tail -n +2  ${results_dir}/UKB_BCT_TWAS_results_${pheno}_bonferroniSigGenes_ref_file >> ${results_dir}/WBC_TWAS_bonferroniSigGenes_nonUniq
	elif [ "$pheno" = "platelet_count" ] || [ "$pheno" = "mean_platelet_volume" ] || [ "$pheno" = "plateletcrit" ] || [ "$pheno" = "platelet_distribution_width" ]; then
        	tail -n +2  ${results_dir}/UKB_BCT_TWAS_results_${pheno}_bonferroniSigGenes_ref_file >> ${results_dir}/PLT_TWAS_bonferroniSigGenes_nonUniq
	else 
		tail -n +2  ${results_dir}/UKB_BCT_TWAS_results_${pheno}_bonferroniSigGenes_ref_file >> ${results_dir}/RBC_TWAS_bonferroniSigGenes_nonUniq
	fi
done


#Determine unique genes across phenotype categories. 
cut -f 1-5 ${results_dir}/WBC_TWAS_bonferroniSigGenes_nonUniq  | sort -k4 | uniq | sed '1 i\chr\tstart_pos\tend_pos\tgene_name\tchunk' > ${results_dir}/WBC_TWAS_bonferroniSigGenes
cut -f 1-5 ${results_dir}/PLT_TWAS_bonferroniSigGenes_nonUniq | sort -k4 | uniq | sed '1 i\chr\tstart_pos\tend_pos\tgene_name\tchunk'> ${results_dir}/PLT_TWAS_bonferroniSigGenes
cut -f 1-5 ${results_dir}/RBC_TWAS_bonferroniSigGenes_nonUniq | sort -k4 | uniq | sed '1 i\chr\tstart_pos\tend_pos\tgene_name\tchunk' > ${results_dir}/RBC_TWAS_bonferroniSigGenes


# #Gene was significant in all studies.
# #awk -v max=$(cut -f 1-5 ${results_dir}/WBC_TWAS_bonferroniSigGenes_nonUniq | sort -k4 | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' | sort -k1nr | head -1 | cut -f 1) '$1==max'  ${results_dir}/WBC_TWAS_bonferroniSigGenes
# 
#Assemble significant gene list, regardless of phenotype category. List of genes which are significant with at least one phenotype - regardless of category.
cat <( head -1 ${results_dir}/WBC_TWAS_bonferroniSigGenes )  <( cat <( tail -n +2 ${results_dir}/WBC_TWAS_bonferroniSigGenes )  <( tail -n +2 ${results_dir}/PLT_TWAS_bonferroniSigGenes ) <( tail -n +2 ${results_dir}/RBC_TWAS_bonferroniSigGenes ) | sort -k4 | uniq ) > ${results_dir}/TWAS_bonferonniSigGenes 

rm ${results_dir}/*_nonUniq
