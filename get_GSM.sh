#!/bin/bash 

## download all of the SRA files associated with a GSM 
GSM=$1

KK=`esearch -db sra -query $GSM | efetch -format runinfo | grep SRR | perl -ne 'm/(SRR\d+).*(SRX\d+)/g; print "$1\t$2\n"'`
SRRS=`echo $KK | perl -ne '@t=split; foreach my $i (@t) {print "$i " if ($i =~ m/SRR/g)}; print "\n"'`
SRXS=`echo $KK | perl -ne '@t=split; foreach my $i (@t) {print "$i " if ($i =~ m/SRX/g)}; print "\n"'`

echo "Getting the following exeriments for $GSM: $KK" 

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
