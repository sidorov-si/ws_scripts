#!/bin/bash 

## can be used for both single and paired-end
## archived FASTQ is assumed

TAG=$1
SPECIES=$2
SINGLE=""
READS=""
REF=""
WDIR=`pwd`

if [[ $TAG == "" || $SPECIES == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag) assuming <tag>.fastq.gz or <tag>.R1.fastq.gz/<tag>.R2.fastq.gz input; 2) species/assembly alias, e.g. genprime_v23"
  exit 1
fi 

if [[ -e $TAG.fastq.gz ]]
then 
  REF="/media/DISK1/reference/kallisto/${SPECIES}_kallisto"
  echo "Processing alignment as single-end, using kallisto index $REF."
  READS="$TAG.fastq.gz"
  SINGLE="--single -l 200 -s 50"
elif [[ -e $TAG.R1.fastq.gz && -e $TAG.R2.fastq.gz ]]
then
  REF="/media/DISK1/reference/kallisto/${SPECIES}_kallisto"
  echo "Processing alignment as paired-end, using kallisto index $REF."
  READS="$TAG.R1.fastq.gz $TAG.R2.fastq.gz"
else
  echo "ERROR: The reqiured fastq.gz files were not found!" 
  exit 1
fi

if [[ ! -e $REF ]]
then 
  echo "ERROR: kallisto index $REF does not exist!" 
  exit 1
fi

kallisto quant -i $REF -t 4 $SINGLE --plaintext -o ${TAG}_kallisto $READS
