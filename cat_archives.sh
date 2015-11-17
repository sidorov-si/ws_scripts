#!/bin/bash

KK=`cat DONOR_IDS`
echo "Listed donors: "$KK
FS=""

for i in $KK
do
  for j in B 12 24
  do
    echo "Processing: $i $j"
    echo "---------------------------------------------"
    for f in ${i}-${j}*.fastq.gz
    do
      FS="$FS ${f%%.gz}"
      echo "inflating $f" 
      gzip -d $f & 
    done
    wait
    echo "now concatenating: $FS"
    cat $FS > ${i}_${j}.fastq
    FS=""
  done
done
