#!/bin/bash 

SPECIES=$1

KK=`for i in *fastq.gz
do 
  TAG1=${i%%.fastq.gz}
  TAG2=${TAG1%%.R?}
  echo $TAG2
done | sort | uniq`

for i in $KK
do 
./kallisto_quant.sh $i $SPECIES &> $i.kallisto_stdout.log 
done
