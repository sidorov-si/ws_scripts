#!/bin/bash 

## download all of the SRA files associated with a GSM 
GSM=$1

echo "Getting the following exeriments for $GSM:" 
esearch -db sra -query $GSM | efetch -format runinfo | awk -F "," '{if (NR>1) printf "%s\t%s\t%s\n",$1,$13,$14}' | grep -v "^\s*$"
SRRS=`esearch -db sra -query $GSM | efetch -format runinfo | awk -F "," '{if (NR>1) print $1}' | grep -v "^\s*$"`
SRXS=`esearch -db sra -query $GSM | efetch -format runinfo | awk -F "," '{if (NR>1) print $13}' | grep -v "^\s*$"`

a=( $SRRS ) 
b=( $SRXS ) 
N=${#a[@]}
M=${#b[@]}

echo "found following srr IDs: "$SRRS", with srx IDs: "$SRXS
for i in `seq 0 $((N-1))`
do
  TRUNC=`echo ${b[$i]} | cut -c 1-6`
  URL="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/$TRUNC/${b[$i]}/${a[$i]}/${a[$i]}.sra"
  echo "downloading from the following URL: $URL"
  wget -b $URL
done
