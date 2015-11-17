#!/bin/bash 

## You need to get 3 files from 
## 1) gene2accession 2) gene2ensembl 3) gene2refseq

## also go to 
## and get "gene_info" files 


zcat gene2accession.gz | grep -P "^9606\t"  > hs_gene2accession & 
zcat gene2accession.gz | grep -P "^10090\t" > mm_gene2accession & 
zcat gene2accession.gz | grep -P "^10116\t" > rn_gene2accession & 
wait

zcat gene2refseq.gz | grep -P "^9606\t"  > hs_gene2refseq & 
zcat gene2refseq.gz | grep -P "^10090\t" > mm_gene2refseq & 
zcat gene2refseq.gz | grep -P "^10116\t" > rn_gene2refseq & 
wait

zcat gene2ensembl.gz | grep -P "^9606\t"  > hs_gene2ensembl & 
zcat gene2ensembl.gz | grep -P "^10090\t" > mm_gene2ensembl & 
zcat gene2ensembl.gz | grep -P "^10116\t" > rn_gene2ensembl & 
wait


