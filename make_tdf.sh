#!/bin/bash

BAM=$1
SPECIES=$2
TAG=${BAM%%.bam}

if [[ $BAM == "" || $SPECIES == "" ]]
then
  echo "ERROR: you must provide both BAM file and species (mm9,mm10,hg19)"
  exit
fi

igvtools count -z 5 -w 50 -e 0 $BAM $TAG.tdf $SPECIES >& $TAG.make_tdf.log 

