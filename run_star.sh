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
./star_align.sh $i $SPECIES &> $i.star_stdout.log 
done
