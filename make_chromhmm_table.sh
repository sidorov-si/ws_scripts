#!/bin/bash

# Systems Immunology Lab, WUSTL

NAME=$1       # common prefix in BED files to be processed, 
              # e.g. [proB]_H3K27ac, [proB]_H3K4me1, etc
WFILE=$2      # genomic windows
CHROMSIZES=$3 # file with chromosome sizes
INDIR=$4      # directory with BED files to be processed
OUTDIR=$5     # output directory for generated tables
TMPDIR=$6     # directory for temp files

echo "This script is intended to run on multiple CPUs"

echo "1. Find intersections between genome windows and marks..."

KK=`awk '{print $1}' $CHROMSIZES`

for i in `ls $INDIR/$NAME*.bed | xargs -I {} sh -c "basename {}"`
do
  bedtools intersect -c -f 0.5 -b $INDIR/$i -a $WFILE > $TMPDIR/${i%%.bed}_binary.bed &
done
wait

echo "2. Generate per-chromosome per-mark binary info..."

for i in `ls $TMPDIR/$NAME*_binary.bed | xargs -I {} sh -c "basename {}"`
do
  for j in $KK
  do
    grep -P "$j\t" $TMPDIR/$i | awk '{print $4}' > $TMPDIR/${i%%_binary.bed}_${j}.x &
  done
  wait
  echo "Processed $i"
done

echo "3. Write out marks..."

ls $TMPDIR/$NAME*_chr1.x | xargs -I {} sh -c "basename {}" | sed "s/_chr1\.x//g" | sed "s/${NAME}_//g" > $TMPDIR/$NAME.$$.marks

echo "Marks are:" `cat $TMPDIR/$NAME.$$.marks`

echo "4. Generate final per-chromosome tables..."

PP=`cat $TMPDIR/$NAME.$$.marks`

for j in $KK
do
  echo "$NAME $j" | sed "s/ /\t/g" > $OUTDIR/${NAME}_${j}_binary.txt
  echo $PP | sed "s/ /\t/g" >>  $OUTDIR/${NAME}_${j}_binary.txt
  paste $TMPDIR/$NAME*_${j}.x >> $OUTDIR/${NAME}_${j}_binary.txt
done

echo "5. Remove temporaty files..."

rm $TMPDIR/{*.x,*binary.bed,$NAME.$$.marks}

echo "Finished."

