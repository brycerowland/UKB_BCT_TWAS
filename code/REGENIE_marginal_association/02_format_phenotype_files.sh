#!/bin/bash

cat <(  head -1 UKB_BCT_adjusted_normalized_phenotypes.tsv | sed 's/ID/FID\tIID/g' ) <( tail -n +2 UKB_BCT_adjusted_normalized_phenotypes.tsv | awk '{ print $1, $0 }' OFS="\t" ) | gzip -c  > UKB_BCT_adjusted_normalized_phenotypes_REGENIE_format.tsv.gz
