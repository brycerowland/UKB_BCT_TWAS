#!/bin/bash

ref_file=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/coloc_assigned_variants_TWAS_misses.tsv
out_dir=/proj/yunligrp/users/bryce/twas/UKB/data/finemapping/coloc_assigned_variants_TWAS_model_liftover

while read chr gene_name; do
	model_dir=/proj/yunligrp/users/bryce/twas/UKB/data/model_training/chr${chr}

	
	for file in ${model_dir}/DGN_chr${chr}_chunk*_${gene_name}_UKB_intersection.betas.EN.txt; do
	if [ ! -f $file ]; then
		continue 
	fi
		cut -f -4 $file | sed 's/^chr//g' | tail -n +2 | awk '{$1=$1} 1' OFS=":" > ${out_dir}/model_b38_ids_chr${chr}_${gene_name}

		~/bin/id_liftOver.sh ${out_dir}/model_b38_ids_chr${chr}_${gene_name} ~/bin/chain_files/hg38ToHg19.over.chain ${out_dir}/liftover_b38_ids_chr${chr}_${gene_name} ${out_dir}/model_b37_ids_chr${chr}_${gene_name}


	done
	
done < <( awk '{split($1,a,":"); print a[1], $3}' $ref_file | tail -n +2 | sort -u  )
