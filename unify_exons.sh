#!/bin/bash 

## make unified coding sequence BED file from GTF file 

GTF=$1


grep -P "\tCDS\t" $GTF  | cut -f 1,4,5 | grep chr | sort -k1,1 -k2,2n > $$.tmp.bed

bedtools merge -i $$.tmp.bed > ${GTF%%.gtf}.merged_CDS.bed

rm $$.tmp.bed
