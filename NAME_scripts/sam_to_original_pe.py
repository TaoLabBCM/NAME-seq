#!/usr/bin/python

import sys
import re
from Bio.Seq import Seq

sam_file = sys.argv[1]
fastq_r1 = sys.argv[2]
fastq_r2 = sys.argv[3]
converted_sam = sam_file
sam_with_original_reads_fwd = sys.argv[4]
sam_with_original_reads_rev = sys.argv[5]
with open(fastq_r1, 'r') as f_r1:
    with open(fastq_r2, 'r') as f_r2:
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
        reads_r2 = {}
        n = 1
        for line in f_r2:
            if n%4 == 1:
                read_id = line.split(' ')[0]
                read_id = read_id[1:]
            if n%4 == 2:
                read_seq = line.rstrip()
                reads_r2[read_id] = read_seq
            n += 1
with open(converted_sam, 'rt') as fin:
    with open(sam_with_original_reads_fwd, 'wt') as fout_fwd:
        with open(sam_with_original_reads_rev, 'wt') as fout_rev:
            n = 1
            for line in fin:
                if line[0] != '@':
                    # print(line)
                    read_id = line.split('\t')[0]
                    # print(read_id)
                    # print(n)'
                    sam_seq = line.split('\t')[9]
                    cigar = line.split('\t')[5]
                    flag = int(line.split('\t')[1])
                    bit_flag = list(bin(flag))[2:]
                    mapq = int(line.split('\t')[4])
                    if len(bit_flag) < 7 or bit_flag[-3] == '1':
                        continue
                    # sam_seq = reads_r1[read_id]
                    if len(bit_flag) >= 7 and bit_flag[-7] == '1':
                        sam_seq = reads_r1[read_id]
                    elif len(bit_flag) >= 7 and bit_flag[-7] == '0':
                        sam_seq = reads_r2[read_id]
                    if len(bit_flag) >= 7 and bit_flag[-5] == '1' and bit_flag[-7] == '1' and bit_flag[-2] == '1' and mapq >= 20:
                        sam_seq = str(Seq(sam_seq).reverse_complement())
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_rev.write(line.replace(line.split('\t')[9], sam_seq))
                    elif len(bit_flag) >= 8 and bit_flag[-5] == '1' and bit_flag[-8] == '1' and bit_flag[-2] == '1' and mapq >= 20:
                        sam_seq = str(Seq(sam_seq).reverse_complement())
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_fwd.write(line.replace(line.split('\t')[9], sam_seq))
                    elif len(bit_flag) >= 7 and bit_flag[-6] == '1' and bit_flag[-7] == '1' and bit_flag[-2] == '1' and mapq >= 20:
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_fwd.write(line.replace(line.split('\t')[9], sam_seq))
                    elif len(bit_flag) >= 8 and bit_flag[-6] == '1' and bit_flag[-8] == '1' and bit_flag[-2] == '1' and mapq >= 20:
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_rev.write(line.replace(line.split('\t')[9], sam_seq))
                else:
                    fout_fwd.write(line)
                    fout_rev.write(line)
                n += 1