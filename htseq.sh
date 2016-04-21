#!/bin/bash

BAM=$1
SPECIES=$2
GTF="/media/DISK1/reference/Annotations/$SPECIES.gtf"
TAG=${BAM%%.bam}


if [[ ! -e $BAM || ! -e $GTF ]]
then
  echo "ERROR: Either $BAM or $GTF file does not exist."
  exit
fi

htseq-count -s no -t exon -f bam $BAM $GTF > $TAG.htseq.out  
