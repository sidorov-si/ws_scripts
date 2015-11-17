#!/bin/bash 

## can be used for both single and paired-end
## archived FASTQ is assumed

TAG=$1
SPECIES=$2
RL=$3
R1=$4
R2=$5
READS=""
WDIR=`pwd`

if [[ $TAG == "" || $SPECIES == "" || $RL == "" || $R1 == "" ]]
then 
  echo "Please provide 1) output name (tag) 2) species 3)read length 4) R1 (or R1/R2 for paired-end) read file name"
  exit 1
fi 

if [[ $R2 == "" ]]
then 
  echo "Processing alignment as single-end, using STAR index /media/DISK1/reference/STAR/${SPECIES}_${RL}bp."
  READS=$WDIR/$R1
else 
  echo "Processing alignment as paired-end, using STAR index /media/DISK1/reference/STAR/${SPECIES}_${RL}bp."
  READS="$WDIR/$R1 $WDIR/$R2"
fi

GENDIR="/media/DISK1/reference/STAR/${SPECIES}_${RL}bp"


mkdir ${TAG}_STAR
cd ${TAG}_STAR
STAR --genomeDir $GENDIR --readFilesIn $READS --runThreadN 4 --readFilesCommand zcat --outFilterMultimapNmax 15 --outFilterMismatchNmax 6  --outSAMstrandField All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM 

mv Aligned.sortedByCoord.out.bam $TAG.bam
mv Aligned.toTranscriptome.out.bam $TAG.tr.bam 
mv Log.out $TAG.star_run.log 
mv Log.final.out $TAG.star_final.log

