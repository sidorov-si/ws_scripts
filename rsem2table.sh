#!/bin/bash 

## this will produce the table suitable for DESeq2 (with integer count values)

NTAG=$1
FILTER=$2
ANN=$3

echo -e "Gene_id\tSymbol\tGene_type" > names
cat $ANN >> names

PP=`ls *genes.results | grep $FILTER`

for i in $PP
do
  echo "processing file $i" 
  TAG=${i%%.genes.results}
  echo $TAG > $TAG.tmp
  #echo $TAG > $TAG.tmp2
  awk '{if (NR>1) print $5}' $i >> $TAG.tmp
  #awk '{if (NR>1) print $7}' $i >> $TAG.tmp2
done

paste names *.tmp  > $NTAG.all_gene.rsem.counts  
#paste names *.tmp2 > $NTAG.all_gene.rsem.fpkms
rm *.tmp names
# rm *.tmp2
