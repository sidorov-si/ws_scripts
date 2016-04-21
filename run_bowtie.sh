#!/bin/bash 

for i in *.fastq.gz
do
  TAG=${i%%.fastq.gz}
  echo $TAG
  ./bowtie2_align.sh $TAG genprime_vM7
done 
