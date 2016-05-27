#!/bin/bash 

SPECIES=$1
STRAND=$2

for i in *tr.bam 
do 
  TAG=${i%%.tr.bam}
  echo "RSEM: processing sample $TAG, file $i"
  ./rsem_quant.sh $TAG $SPECIES $STRAND
done
