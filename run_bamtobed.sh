#!/bin/bash

BAM=$1
TAG=${BAM%%.bam}

bedtools bamtobed -i $BAM > $TAG.bed 
