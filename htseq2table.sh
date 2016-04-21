#!/bin/bash 

## this will produce the table suitable for DESeq2 (with integer count values)

NTAG=$1
FILTER=$2
ANN=$3

echo "using annotation $ANN"

PP=`ls *htseq.out | grep $FILTER`

for i in $PP
do
  echo "processing file $i" 
  TAG=${i%%.htseq.out}
  echo $TAG > $TAG.tmp
  grep -v -P "__" $i | awk '{print $2}' >> $TAG.tmp
done

paste $ANN *.tmp > $NTAG.all_gene.htseq.counts
rm *.tmp 
