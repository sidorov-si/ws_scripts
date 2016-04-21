#!/bin/bash 

## can be used for both single and paired-end
## archived FASTQ is assumed

TAG=$1
SPECIES=$2
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
  RL="101"
  REF="/media/DISK1/reference/STAR/${SPECIES}_${RL}bp"
  echo "Processing alignment as single-end, using STAR index $REF."
  READS=$WDIR/$TAG.fastq.gz
else
  echo "ERROR: The reqiured fastq.gz files were not found!" 
  exit 1
fi

if [[ ! -d $REF ]]
then 
  echo "ERROR: STAR index $REF does not exist!" 
  exit 1
fi

mkdir ${TAG}_STAR
cd ${TAG}_STAR
STAR --genomeDir $REF --readFilesIn $READS --runThreadN 4 --readFilesCommand zcat --outFilterMultimapNmax 15 --outReadsUnmapped Fastx --seedSearchStartLmax 30 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 50  --outSAMstrandField All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM 

mv Aligned.sortedByCoord.out.bam $TAG.bam
mv Aligned.toTranscriptome.out.bam $TAG.tr.bam 
mv Log.out $TAG.star_run.log 
mv Log.final.out $TAG.star_final.log

