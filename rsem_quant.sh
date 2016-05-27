#!/bin/bash 

## can be used for both single and paired-end
## STAR transcriptome BAM is assumed to be named $TAG.tr.bam 


TAG=$1
SPECIES=$2
STRAND=$3
PARROT=""
TEST=""
REF=/media/DISK1/reference/RSEM/${SPECIES}_rsem
WDIR=`pwd`

if [[ $TAG == "" || $SPECIES == "" || $STRAND == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag) assuming <tag>.fastq.gz or <tag>.R1.fastq.gz/<tag>.R2.fastq.gz input; 2) species/assembly alias, e.g. genprime_v23; 3) strandedness as NONE/FR/RF"
  exit 1
fi 

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

TEST=`samtools view -f 0x1 $TAG.tr.bam | head`

if [[ $TEST == "" ]]
then
  PARROT="" 
else
  PARROT="--paired-end"
fi

if [[ ! -e $REF.grp ]]
then 
  echo "ERROR: rsem index $REF does not exist!" 
  exit 1
fi

echo "Executing command: rsem-calculate-expression -p 4 --bam --no-bam-output $HTSSTR --estimate-rspd $PARROT $TAG.tr.bam $REF $TAG &> $TAG.rsem.log"
rsem-calculate-expression -p 4 --bam --no-bam-output $HTSSTR --estimate-rspd $PARROT $TAG.tr.bam $REF $TAG &> $TAG.rsem.log 
