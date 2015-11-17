#!/bin/bash

for i in *sra
do
  fastq-dump --split-3 $i & 
done

wait

for i in *fastq
do
  gzip $i &
done

wait
