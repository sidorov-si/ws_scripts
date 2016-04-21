#!/bin/bash 

## this will produce the table suitable for DESeq2 (with integer count values)

NTAG=$1
FILTER=$2

F1=`ls *genes.results | head -n 1`
awk '{print $1}' $F1 > names 
PP=`ls *genes.results | grep $FILTER`

for i in $PP
do
  echo "processing file $i" 
  TAG=${i%%.genes.results}
  echo $TAG > $TAG.tmp
  echo $TAG > $TAG.tmp2
  awk '{if (NR>1) print $5}' $i >> $TAG.tmp
done

paste names *.tmp  > $NTAG.all_gene.rsem.counts  
rm *.tmp names
