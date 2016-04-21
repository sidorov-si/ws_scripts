#!/bin/bash 

## tab-delimited list of GSM-name correspondences (no spaces!)
SMP=$1

KK=`cat $SMP`
N=`wc -l $SMP | awk '{print $1}'`
TYPE="undef"

a=( $KK ) 

for i in `seq 0 2 $((2*N-2))`
do
  j=$((i+1))
  ##echo "$i $j"
  SRR=`esearch -db sra -query "${a[$i]}" | efetch -format runinfo | grep SRR | awk -F "," '{print $1}'`
  b=( $SRR ) 
  if [[ -e ${b[0]}_1.fastq.gz ]]
  then 
    TYPE="paired"
  elif [[ -e ${b[0]}.fastq.gz ]]
  then
    TYPE="single"
  else 
    TYPE="undef"
  fi 
  
  echo "Processing file: records: ${a[$i]} ${a[$j]}, testing file: ${b[0]}, experiment type: $TYPE, SRR list: $SRR"
  if [[ $TYPE == "single" ]]
  then
    SRR_FQ=`echo $SRR | perl -ne 's/SRR(\d+)/SRR$1.fastq.gz/g; print'`
    
    if [[ "$SRR_FQ" = "${SRR_FQ%[[:space:]]*}" ]] 
    then
      echo "will be MOVING single file  $SRR_FQ into ${a[$j]}.fastq.gz"
      mv $SRR_FQ  ${a[$j]}.fastq.gz
    else 
      echo "will be CONCATENATING files $SRR_FQ into ${a[$j]}.fastq.gz"
      zcat $SRR_FQ | gzip > ${a[$j]}.fastq.gz
    fi

  elif [[ $TYPE == "paired" ]]
  then
    SRR_FQ1=`echo $SRR | perl -ne 's/SRR(\d+)/SRR$1_1.fastq.gz/g; print'`
    SRR_FQ2=`echo $SRR | perl -ne 's/SRR(\d+)/SRR$1_2.fastq.gz/g; print'`
  
    if [[ "$SRR_FQ1" = "${SRR_FQ1%[[:space:]]*}" ]] 
    then
      echo "will be MOVING single file $SRR_FQ1 to ${a[$j]}.R1.fastq.gz"
      echo "will be MOVING single file $SRR_FQ2 to ${a[$j]}.R2.fastq.gz"
      mv $SRR_FQ1 ${a[$j]}.R1.fastq.gz
      mv $SRR_FQ2 ${a[$j]}.R2.fastq.gz
    else 
      echo "will be CONCATENATING files $SRR_FQ1 to ${a[$j]}.R1.fastq"
      echo "will be CONCATENATING files $SRR_FQ2 to ${a[$j]}.R2.fastq"
      zcat $SRR_FQ1 | gzip > ${a[$j]}.R1.fastq.gz
      zcat $SRR_FQ2 | gzip > ${a[$j]}.R2.fastq.gz
    fi
  else 
    echo "ERROR: something went terribly wrong! Your experiment is neither paired nor single!"
  fi
done   

