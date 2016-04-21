#!/bin/bash 

## http or ftp - makes virtually no difference.

KK=`cat SAMPLES`

for i in $KK
do
  URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/$i"
  echo $i
  wget $URL 
done 

