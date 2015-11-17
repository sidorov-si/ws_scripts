#!/bin/bash 


SAMPLES=$1  ## file with sample dirs 
KK=`cat $SAMPLES`

for i in $KK
do
  echo "downloading files from folder $i"
  wget -b -r -nd --no-parent http://sysg1.cs.yale.edu:3010/1oW9aOLkIK0d5PaQjCVfcRQn743Lc/Project_Os83/$i
  echo "now sleeping for 10 min" 
  sleep 600
done
