#!/bin/bash 


TAG=$1
SPECIES=$2 ## use smth like genprime_v23 etc

if [[ $TAG == "" || $SPECIES == "" ]]
then 
  echo "Please provide 1) output name (tag) 2) species alias 3) read R1 (or R1/R2 for paired-end) file name"
  exit 1
fi 

if [[ -e $TAG.fastq.gz ]]
then
  echo "Aligning single-end reads ($TAG.fastq.gz), --sensitive-local mode, against index $SPECIES, tag $TAG"
  bowtie2 --sensitive-local -t -p 4 -S $TAG.sam -x /media/DISK1/reference/Bowtie2/$SPECIES -U $TAG.fastq.gz &> $TAG.bowtie2.log  
  samtools view -bS -q 10 $TAG.sam > $TAG.bam
  samtools sort -@ 4 -T $TAG -o $TAG.sorted.bam $TAG.bam
  mv $TAG.sorted.bam $TAG.bam 
  rm $TAG.sam
elif [[ -e $TAG.R1.fastq.gz && -e $TAG.R2.fastq.gz ]]
then
  echo "Aligning paired-end reads ($TAG.R1.fastq.gz,$TAG.R2.fastq.gz), --sensitive-local mode, against index $SPECIES, tag $TAG"
  bowtie2 --sensitive-local -t -p 4 -S $TAG.sam -x /media/DISK1/reference/Bowtie2/$SPECIES -1 $TAG.R1.fastq.gz -2 $TAG.R2.fastq.gz &> $TAG.bowtie2.log
  samtools view -bS -q 10 $TAG.sam > $TAG.bam
  samtools sort -@ 4 -T $TAG -o $TAG.sorted.bam $TAG.bam
  mv $TAG.sorted.bam $TAG.bam 
  rm $TAG.sam
else 
  echo "The tag $TAG cannot be classified as either single- or paired-end!"
fi
