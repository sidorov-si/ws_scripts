#!/bin/bash

LIST=$1
N=`wc -l $LIST | awk '{print $1-1}'`
KK=`awk '{print $1}' $LIST`
PP=`awk '{print $2}' $LIST`
a=( $KK )
b=( $PP ) 

for i in `seq 0 $N`
do
  mkdir ${a[$i]}_SICER
  cd ${a[$i]}_SICER
  echo "Processing TAG ${a[$i]} with gap size ${b[$i]} bp" 
  nohup SICER-rb.sh /media/DISK2/ssidorov/mESC/selected ${a[$i]}.bed . mm10 1 200 250 0.77 ${b[$i]} 100  &> ${a[$i]}.sicer.log &  
  cd ..
done
