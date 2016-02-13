#!/usr/bin/env python

import sys

def find_name(gene_id, all_genes_lines):
    for line in all_genes_lines:
        if gene_id in line:
            gene_name = line.split('\t')[12]
    return gene_name

# GTF file with genes and their IDs (both gene_id and transcription_id)
gtf_in_filename = sys.argv[1] 
# File with genes and all info about them including IDs and names
all_genes_filename = sys.argv[2] 
# Output file with a list of gene names.
# In case several loci have the same gene name, their names are
# augmented by coordinates
out_filename = '.'.join(gtf_in_filename.split('.')[:-1]) + '.names.txt' 

with open(all_genes_filename, 'r') as all_genes_file:
    all_genes_lines = all_genes_file.readlines()

with open(gtf_in_filename, 'r') as infile, open(out_filename, 'w') as outfile:
    gene_name = ''
    gene_id = ''
    gene_start = ''
    gene_end = ''
    gene_chr = ''
    print 'Begin...'
    for num, line in enumerate(infile):
        if num % 100 == 0 and num != 0:
            print 'Processed ' + str(num) + ' genes.'
        gtf_line = line.strip().split('\t')
        gene_id = (gtf_line[8].split('"'))[1]
        gene_start = gtf_line[3]
        gene_end = gtf_line[4]
        gene_chr = gtf_line[0]
        gene_name = find_name(gene_id, all_genes_lines)
        outfile.write(gene_name + '\t' + gene_chr + ':' + gene_start + '-' + gene_end + \
                      '\t' + gene_id + '\n')

    print 'End.'

