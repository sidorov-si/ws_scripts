#!/bin/bash

for tag in `ls *.fastq.gz | awk -F"." '{print $1}' | tr '\n' ' '` 
do 
    ./bowtie2_align.sh $tag hg19 $tag.fastq.gz
done

