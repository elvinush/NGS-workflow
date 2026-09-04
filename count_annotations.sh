#!/bin/bash
# counts how many ANN entries snpEff put in a vcf. one variant can have several
# annotations so this is not the same as the number of variants.
# usage: ./count_annotations.sh annotation/3GC_genome_annotation.vcf.gz

VCF=$1

zgrep -v "^#" $VCF \
| awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^ANN=/) print $i}' \
| sed 's/^ANN=//' \
| tr ',' '\n' \
| wc -l
