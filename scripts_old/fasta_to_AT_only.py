#!/usr/bin/python

import sys
import re
import os

fasta_file = sys.argv[1]
fasta_file_AT_only = sys.argv[2]

if not os.path.exists(os.path.dirname(fasta_file_AT_only)):
    try:
        os.makedirs(os.path.dirname(fasta_file_AT_only))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
            
with open(fasta_file,'r') as fasta:
    with open(fasta_file_AT_only,'w') as fasta_out:
        for line in fasta:
            if line[0] == '>':
                fasta_out.write(line)
            if line[0] != '>':
                line = str(line).upper()
                line = line.replace('G','A')
                line = line.replace('C','T')
                fasta_out.write(line)