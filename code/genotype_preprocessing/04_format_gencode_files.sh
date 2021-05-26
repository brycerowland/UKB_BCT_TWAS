#!/bin/bash

#Performs training by chromsome
loop_chr=$1

data_dir="/proj/yunligrp/users/bryce/twas/UKB/data"
out_dir=${data_dir}/predicted_gexp/gene_ref_files/chr${loop_chr}

rm -f chunk_vector_chr${loop_chr}

#Find the corresponding chunk for each gene.
while read chr start_pos end_pos gene_name
do
	window_start=$(($start_pos-1000000<1?1:$start_pos-1000000))
        window_end=$(( end_pos + 1000000 ))

        #Find a chunk that works.
        for i in {1..42};
        do
                if [ "$window_start" -ge $(( ($i-1)*6000000 + 1)) ] && [ "$window_start" -lt $(( $i*6000000 + 4000000 )) ] && [ "$window_end" -ge $(( ($i-1)*6000000 + 1)) ] && [ "$window_end" -lt $(( $i*6000000 + 4000000 )) ]
                then
                        chunk=$i
                fi
        done

        #Update reference file to include DGN_chunk
        echo $chunk >> chunk_vector_chr${chr}
	

done < <( awk -v chr=$loop_chr '$1=="chr"chr' ${data_dir}/DGN_subset_gencode_build38_genes | sed 's/^chr//g' )

paste <( awk -v chr=$loop_chr '$1=="chr"chr' ${data_dir}/DGN_subset_gencode_build38_genes| sed 's/^chr//g' ) chunk_vector_chr${loop_chr} > ${out_dir}/gencode_build38_chr${loop_chr}_genes.txt
