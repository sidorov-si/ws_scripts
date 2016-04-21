#!/bin/bash

## no gzipping, large datasets, 10 fq-dumps a turn

ls *sra > LIST
split -l 10 -d LIST

for i in x??
do
  KK=`cat $i`
  for j in $KK
  do
    fastq-dump --split-3 $j & 
  done
  wait
done
