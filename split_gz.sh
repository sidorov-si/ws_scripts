#!/bin/bash

## no gzipping, large datasets, 10 fq-dumps a turn

ls *fastq > LIST
split -l 8 -d LIST

for i in x??
do
  KK=`cat $i`
  for j in $KK
  do
    gzip $j & 
  done
  wait
done
