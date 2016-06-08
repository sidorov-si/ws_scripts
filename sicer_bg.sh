#!/bin/bash

## one-sample script to run SICER.sh with Input/control
## run via run_sicer.sh 

BEDDIR=$1
BED=$2
INPBED=$3
SPECIES=$4
  
TAG=${BED%%.bed}
mkdir ${TAG}_SICER
cd ${TAG}_SICER
GAP=200

if [[ $TAG == *H3K9me* || $TAG == *H3K27me* || $TAG == *H3K36me* || $TAG == *H3K79me* ]]
then
  GAP=600
fi

echo "SICER.sh: processing TAG $TAG with gap size $GAP bp and input file $INPBED." 
SICER.sh $BEDDIR $BED $INPBED . $SPECIES 1 200 250 0.77 $GAP 100  &> $TAG.sicer.log 
