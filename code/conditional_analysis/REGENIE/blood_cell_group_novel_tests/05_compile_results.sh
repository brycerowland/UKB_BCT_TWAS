#!/bin/bash

data_dir=/proj/yunligrp/users/bryce/twas/UKB/data/conditional_analysis/REGENIE/blood_cell_group_novel_tests/regenie_step2


for file in ${data_dir}/*regenie; do

	gene_name=$( basename $file | sed 's/UKB_TWAS_conditional_test_novelByClass_//g' | cut -f 1 -d '_' )

	grep $gene_name $file

done > results.txt
