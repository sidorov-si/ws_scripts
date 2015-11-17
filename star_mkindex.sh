#!/bin/bash 


## this is a simplified 1-step generation command
## WARNING: you still have to check if the junctions file have junctions that are < 20 bp or esp. of negative length. 

FA=$1
GTF=$2
TAG=$3
RL=$4 ## read length

if [[ $FA == "" || $GTF == "" || $TAG == "" || $RL = "" ]]
then 
  echo "Please provide reference genome sequence, annotation, output name (tag), and read length"
  exit
fi 

mkdir ${TAG}_${RL}bp
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir ${TAG}_${RL}bp --genomeFastaFiles $FA --sjdbGTFfile $GTF --sjdbOverhang $((RL-1)) >& ${TAG}_${RL}bp.star_mkindex.log  


