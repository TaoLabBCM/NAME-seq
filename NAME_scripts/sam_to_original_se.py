#!/usr/bin/python

import sys
import re
from Bio.Seq import Seq

sam_file = sys.argv[1]
fastq_r1 = sys.argv[2]

converted_sam = sam_file
sam_with_original_reads_fwd = sys.argv[3]
sam_with_original_reads_rev = sys.argv[4]

with open(fastq_r1, 'r') as f_r1:
    reads_r1 = {}
    n = 1
    for line in f_r1:
        if n%4 == 1:
            read_id = line.split(' ')[0]
            read_id = read_id[1:]
        if n%4 == 2:
            read_seq = line.rstrip()
            reads_r1[read_id] = read_seq
        n += 1

with open(converted_sam, 'rt') as fin:
    with open(sam_with_original_reads_fwd, 'wt') as fout_fwd:
        with open(sam_with_original_reads_rev, 'wt') as fout_rev:
            n = 1
            for line in fin:
                if line[0] != '@':
                    line_list = line.split('\t')
                    read_id = line_list[0]
                    sam_seq = line_list[9]
                    cigar = line_list[5]
                    mapq = int(line_list[4])
                    flag = int(line_list[1])
                    if flag == 16 and mapq >= 20:
                        sam_seq = reads_r1[read_id]
                        sam_seq = str(Seq(sam_seq).reverse_complement())
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_rev.write(line.replace(line_list[9], sam_seq))
                    elif flag == 0 and mapq >= 20:
                        sam_seq = reads_r1[read_id]
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_fwd.write(line.replace(line_list[9], sam_seq))
                else:
                    fout_fwd.write(line)
                    fout_rev.write(line)