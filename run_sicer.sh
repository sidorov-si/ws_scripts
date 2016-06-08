#!/bin/bash 

## put this and sicer_nobg.sh into input BED directory, and run it here. 
BEDDIR=`pwd`
SPECIES=$1

for i in *bed
do 
  ./sicer_nobg.sh $BEDDIR $i $SPECIES 
done

for i in *SICER
do 
  
