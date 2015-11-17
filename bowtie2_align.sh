#!/bin/bash 


TAG=$1
SPECIES=$2 ## use smth like genprime_v23 etc
R1=$3  ## archives are OK 
R2=$4

if [[ $TAG == "" || $SPECIES == "" || $R1 == "" ]]
then 
  echo "Please provide 1) output name (tag) 2) species alias 3) read R1 (or R1/R2 for paired-end) file name"
  exit 1
fi 

if [[ $R4 == "" ]]
then
  echo "Aligning single-end reads ($R1), --sensitive-local mode, against index $SPECIES, tag $TAG"
  bowtie2 --sensitive-local -t -p 4 -S $TAG.sam -x /media/DISK1/reference/Bowtie2/$SPECIES -U $R1 
  samtools view -bS -q 10 $TAG.sam > $TAG.bam
  samtools sort $TAG.bam $TAG.sorted
  mv $TAG.sorted.bam $TAG.bam 
  rm $TAG.sam
else 
  echo "Aligning paired-end reads ($R1,$R2), --sensitive-local mode, against index $SPECIES, tag $TAG"
  bowtie2 --sensitive-local -t -p 4 -S $TAG.sam -x /media/DISK1/reference/Bowtie2/$SPECIES -1 $R1 -2 $R2
  samtools view -bS -q 10 $TAG.sam > $TAG.bam
  samtools sort $TAG.bam $TAG.sorted
  mv $TAG.sorted.bam $TAG.bam 
  rm $TAG.sam
fi
