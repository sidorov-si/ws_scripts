#!/bin/bash 

RMSK=$1
GENCODE=$2
HEADER=$3
TAG=$4 ## eg gencode_vM6 

grep rRNA $RMSK | awk -F "\t" '{if ($6 !~ "_") print}' > rRNA_from_rmsk.table  ## we don't want any patches/scaffolds here since they won't be in Gencode notation. 
grep -P "\tgene\t" $GENCODE | grep rRNA > rRNA_from_gencode.gtf 

awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\n",$6,$7,$8,$11,$12,$10}' rRNA_from_rmsk.table | sort -k1,1V -k2,2n -k3,3n > rRNA_from_rmsk.bed
perl -ne 'm/(chr.*?)\tENSEMBL\tgene\t(\d+?)\t(\d+?)\t\.\t(.*?)\t.*gene_id \"(.*?)\".*gene_name \"(.*?)\"/; print "$1\t$2\t$3\t$6\t$5\t$4\n"' rRNA_from_gencode.gtf | sort -k1,1V -k2,2n -k3,3n > rRNA_from_gencode.bed

KK1=`awk '{sum+=$3-$2} END {print sum}' rRNA_from_rmsk.bed`
KK2=`awk '{sum+=$3-$2} END {print sum}' rRNA_from_gencode.bed`
KK3=`bedtools intersect -wao -a rRNA_from_gencode.bed -b rRNA_from_rmsk.bed | awk '{sum+=$13} END {print sum}'`

echo "overall coverage: $KK1 bp in rmsk intervals, $KK2 bp in gencode inervals, $KK3 bp common"

cat rRNA_from_gencode.bed rRNA_from_rmsk.bed | sort -k1,1V -k2,2n -k3,3n > rRNA_combined.bed
bedtools merge -c 4 -o collapse -s -i rRNA_combined.bed > rRNA_merged.bed
KK4=`awk '{sum+=$3-$2} END {print sum}' rRNA_merged.bed`

echo "coverage of the merged intervals file: $KK4 bp"
cat $HEADER rRNA_merged.bed > $TAG.rRNA_merged.intervals
