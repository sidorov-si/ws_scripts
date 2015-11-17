#!/bin/bash

## AP 

FA=$1
TAG=$2

if [[ $FA == "" || $TAG == "" ]]
then
  echo "You must provide both reference FASTA and the base nametag (ie mm10 or gencode_vM6)"
  exit
fi

bowtie-build $FA $TAG >& $TAG.bowtie_build.log 
