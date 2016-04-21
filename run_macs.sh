#!/bin/bash 

SPECIES=$1

if [[ $SPECIES != "mm" && $SPECIES != "hs" ]]
then
  echo "ERROR: Species must be one of the following: hs, mm, ce, dm, or just use a number for effective genome size." 
  exit 1
fi 

for i in *bam
do
  echo "Processing file $i with NO BACKGROUND experiment."
  TAG=${i%%.bam}
  nohup macs14 -t $i -f BAM -g hs -n $TAG --verbose=2  &> $TAG.macs.log & 
done 
