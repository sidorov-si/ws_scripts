#!/bin/bash

SPECIES=$1
STRAND=$2
HTSSTR=""
REF=""

if [[ $STRAND == "NONE" ]]
then
  HTSSTR="--forward-prob 0.5"
  echo "Proceeding using stranded setting $STRAND"
elif [[ $STRAND == "FR" ]] 
then
  HTSSTR="--forward-prob 1.0"
  echo "Proceeding using stranded setting $STRAND"
elif [[ $STRAND == "RF" ]]
then
  HTSSTR="--forward-prob 0.0"
  echo "Proceeding using stranded setting $STRAND"
else 
  echo "ERROR: you must set strand variable to either NONE, FR, or RF"
  exit
fi

REF="/media/DISK1/reference/RSEM/${SPECIES}_rsem"

for i in *tr.bam
do
  TAG=${i%%.tr.bam}
  echo "Processing file $i, output name $TAG"
  rsem-calculate-expression -p 4 --bam --no-bam-output $HTSSTR --estimate-rspd $i $REF $TAG &> $TAG.rsem.log 
done
