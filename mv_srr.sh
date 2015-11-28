#!/bin/bash 

## add a mark to SRR file names corresponding to the given GSM

mark=$1
GSM=$2
ext=$3

echo "Getting the following exeriments for $GSM:" 
esearch -db sra -query $GSM | efetch -format runinfo | awk -F "," '{if (NR>1) printf "%s\t%s\t%s\n",$1,$13,$14}' | grep -v "^\s*$"
SRRS=`esearch -db sra -query $GSM | efetch -format runinfo | awk -F "," '{if (NR>1) print $1}' | grep -v "^\s*$"`
SRXS=`esearch -db sra -query $GSM | efetch -format runinfo | awk -F "," '{if (NR>1) print $13}' | grep -v "^\s*$"`

a=( $SRRS ) 
b=( $SRXS ) 
N=${#a[@]}
M=${#b[@]}

for i in `seq 0 $((N-1))`
do
     mv "${a[$i]}.$ext" "$mark.${a[$i]}.$ext"
done   

