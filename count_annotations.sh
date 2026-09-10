#!/bin/bash
# counts how many ANN entries snpEff put in a vcf. one variant can have several
# annotations so this is not the same as the number of variants.
# usage: ./count_annotations.sh annotation/3GC_genome_annotation.vcf.gz

VCF=$1

zgrep -v "^#" $VCF \
| cut -f8 \
| grep -o 'ANN=[^;]*' \
| sed 's/^ANN=//' \
| tr ',' '\n' \
| wc -l
