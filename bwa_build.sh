#!/bin/bash

## AP 

FA=$1
TAG=$2

if [[ $FA == "" || $TAG == "" ]]
then
  echo "You must provide both reference FASTA and the base nametag (ie mm10 or gencode_vM6)"
  exit
fi

WDIR=`pwd`
mkdir ${TAG}_bwa_index
cd ${TAG}_bwa_index
bwa index $WDIR/$FA  >& $TAG.bwa_index.log 
