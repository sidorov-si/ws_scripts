#!/usr/bin/env python

import sys
import subprocess
# GTF file with start- and stop-codons, exons, and CDSs
gtf_in_filename = sys.argv[1] 
# GTF files with whole genes as features
gtf_out_filename = '.'.join(gtf_in_filename.split('.')[:-1]) + '.genes.gtf' 
with open(gtf_in_filename, 'r') as ingtf, open(gtf_out_filename, 'w') as outgtf:
    gene_id = ''
    gene_start = 0
    gene_end = 0
    print 'Begin...'
    for num, line in enumerate(ingtf):
        if num % 1000 == 0 and num != 0:
            print 'Processed ' + str(num) + ' features.'
        gtf_line = line.strip().split('\t')
        gtf_feature = gtf_line[2]
        if gtf_feature == 'exon':
            current_gene_id = (gtf_line[8].split('"'))[1]
            if gene_id != current_gene_id:
                if gene_end != 0:
                    outgtf.write(chr_name + '\t' + ref_seq + '\t' + feature + '\t' + \
                                 gene_start + '\t' + gene_end + '\t' + score + '\t' + \
                                 strand + '\t' + codon_position + '\t' + \
                                 'gene_id ' + '"' + gene_id  + '"\n')
                gene_id = current_gene_id
                chr_name = gtf_line[0]
                ref_seq = gtf_line[1]
                feature = 'gene'
                gene_start = gtf_line[3]
                gene_end = gtf_line[4]
                score = gtf_line[5]
                strand = gtf_line[6]
                codon_position = gtf_line[7]
            else:
                gene_end = gtf_line[4]
                
    outgtf.write(chr_name + '\t' + ref_seq + '\t' + feature + '\t' + \
                 gene_start + '\t' + gene_end + '\t' + score + '\t' + \
                 strand + '\t' + codon_position + '\t' + \
                 'gene_id ' + '"' + gene_id  + '"\n')
    print 'End.'

