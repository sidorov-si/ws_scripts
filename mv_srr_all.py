#!/usr/bin/env python

import sys
import subprocess

gsm_mark_list = sys.argv[1]
ext = sys.argv[2]
with open(gsm_mark_list, 'r') as infile:
    for line in infile:
        mark, gsm = line.strip().split('\t')
        subprocess.Popen(['./mv_srr.sh', mark, gsm, ext])

