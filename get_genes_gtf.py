#!/usr/bin/env python

import sys
# GTF file with exons;
# this file should be sorted by (in the following order!): 
# gene_id, transcript_id, chromosome name, strand, and start coordinate 
gtf_in_filename = sys.argv[1] 
# GTF file with whole genes as features
gtf_out_filename = '.'.join(gtf_in_filename.split('.')[:-1]) + '.genes.gtf' 
with open(gtf_in_filename, 'r') as ingtf, open(gtf_out_filename, 'w') as outgtf:
    gene_id = ''
    transcript_id = ''
    gene_start = -1
    gene_end = -1
    chr_name = ''
    strand = ''
    print 'Begin...'
    for num, line in enumerate(ingtf):
        if num % 1000 == 0 and num != 0:
            print 'Processed ' + str(num) + ' features.'
        gtf_line = line.strip().split('\t')
        gtf_feature = gtf_line[2]
        current_gene_id = (gtf_line[8].split('"'))[1]
        current_transcript_id = (gtf_line[8].split('"'))[3]
        current_chr_name = gtf_line[0]
        current_strand = gtf_line[6]
        if gene_id != current_gene_id or transcript_id != current_transcript_id or \
           chr_name != current_chr_name or strand != current_strand:
            if gene_end != -1:
                outgtf.write(chr_name + '\t' + ref_seq + '\t' + 'gene' + '\t' + \
                             str(gene_start) + '\t' + str(gene_end) + '\t' + score + '\t' + \
                             strand + '\t' + codon_position + '\t' + \
                             'gene_id ' + '"' + gene_id  + '" ' \
                             'transcript_id ' + '"' + transcript_id  + '"\n')
            gene_id = current_gene_id
            transcript_id = current_transcript_id
            chr_name = gtf_line[0]
            ref_seq = gtf_line[1]
            gene_start = min(int(gtf_line[3]), int(gtf_line[4]))
            gene_end = max(int(gtf_line[3]), int(gtf_line[4]))
            score = gtf_line[5]
            strand = gtf_line[6]
            codon_position = gtf_line[7]                
        else:
            min_coord = min(int(gtf_line[3]), int(gtf_line[4]))
            if min_coord < gene_start:
                gene_start = min_coord
            max_coord = max(int(gtf_line[3]), int(gtf_line[4]))
            if max_coord > gene_end:
                gene_end = max_coord

    outgtf.write(chr_name + '\t' + ref_seq + '\t' + 'gene' + '\t' + \
                 str(gene_start) + '\t' + str(gene_end) + '\t' + score + '\t' + \
                 strand + '\t' + codon_position + '\t' + \
                 'gene_id ' + '"' + gene_id  + '" ' \
                 'transcript_id ' + '"' + transcript_id  + '"\n')
               
    print 'End.'

