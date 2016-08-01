#!/usr/bin/env bash
# Genes from Gencode.
genes=$1
# Output file suitable for Picard CollectRnaSeqMetrics.jar.
rRNA=$2

# Intervals for rRNA transcripts.
grep -E 'gene_type "rRNA"|gene_type "Mt_rRNA"' $genes | \
    awk '$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> $rRNA

