#!/usr/bin/env python

import sys
# GTF file with exons;
# this file should be sorted by (in the following order!): 
# gene_id, transcript_id, chromosome name, strand, and start coordinate 
gtf_in_filename = sys.argv[1] 
# GTF file with introns as features
gtf_out_filename = '.'.join(gtf_in_filename.split('.')[:-1]) + '.introns.gtf' 
with open(gtf_in_filename, 'r') as ingtf, open(gtf_out_filename, 'w') as outgtf:
    gene_id = ''
    transcript_id = ''
    chr_name = ''
    strand = ''
    prev_exon_start = -1
    prev_exon_end = -1
    print 'Begin...'
    for num, line in enumerate(ingtf):
        if num % 1000 == 0 and num != 0:
            print 'Processed ' + str(num) + ' features.'
        gtf_line = line.strip().split('\t')
        current_gene_id = (gtf_line[8].split('"'))[1]
        current_transcript_id = (gtf_line[8].split('"'))[3]
        current_chr_name = gtf_line[0]
        current_strand = gtf_line[6]
        current_exon_start = min(int(gtf_line[3]), int(gtf_line[4]))
        current_exon_end = max(int(gtf_line[3]), int(gtf_line[4]))
        if gene_id != current_gene_id or transcript_id != current_transcript_id or \
           chr_name != current_chr_name or strand != current_strand:
            prev_exon_end = current_exon_end
            gene_id = current_gene_id
            transcript_id = current_transcript_id
            chr_name = current_chr_name
            strand = current_strand
            continue
        intron_start = prev_exon_end + 1
        intron_end = current_exon_start - 1
        prev_exon_end = current_exon_end
        source = gtf_line[1]
        outgtf.write(chr_name + '\t' + source + '\t' + 'intron' + '\t' + \
                     str(intron_start) + '\t' + str(intron_end) + '\t' + '0.0' + '\t' + \
                     strand + '\t' + '.' + '\t' + 'gene_id ' + '"' + gene_id  + '" ' \
                     'transcript_id ' + '"' + transcript_id  + '"\n')
        gene_id = current_gene_id
        transcript_id = current_transcript_id
        chr_name = current_chr_name
        strand = current_strand

    print 'End.'

