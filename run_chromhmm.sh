#/bin/bash

for i in hESC_input  HUVEC_input  IMR90_input
do
  for j in `seq 26 40`
  do
    TAG=${i%%_input}
    echo "Processing cell type $TAG, with $j states"
    java -mx30G -jar ~/ChromHMM/ChromHMM.jar LearnModel -p 4 $i ${TAG}_${j}.output $j hg19 > ${TAG}_${j}_chromhmm.log & 
    wait 
  done
done
