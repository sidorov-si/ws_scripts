#!/bin/bash 

GSM_LIST=$1

for i in `cat $GSM_LIST`
do
  ./get_GSM.sh $i
  sleep 180
done 
