#!/usr/bin/python

import sys
import re
from Bio.Seq import Seq
from Bio import SeqIO

# get sam file, fastq files, output files
sam_file = sys.argv[1]
fastq_r1 = sys.argv[2]
fastq_r2 = sys.argv[3]
converted_sam = sam_file
sam_with_original_reads_fwd = sys.argv[4]
sam_with_original_reads_rev = sys.argv[5]

# loading all fastq to dict takes too much memory. Instead, we can create index file for fastq using SeqIO.index to query each read.
index_r1 = SeqIO.index(fastq_r1, "fastq")  # Use "fasta" if your file is in FASTA format
index_r2 = SeqIO.index(fastq_r2, "fastq")  # Use "fasta" if your file is in FASTA format

# converted sam to original sam
with open(converted_sam, 'rt') as fin:
    with open(sam_with_original_reads_fwd, 'wt') as fout_fwd:
        with open(sam_with_original_reads_rev, 'wt') as fout_rev:
            for line in fin:
                if line[0] != '@':
                    # load read id, sequence, cigar, flag, and mapq from sam files
                    read_id = line.split('\t')[0]
                    sam_seq = line.split('\t')[9]
                    cigar = line.split('\t')[5]
                    flag = int(line.split('\t')[1])
                    bit_flag = list(bin(flag))[2:]
                    mapq = int(line.split('\t')[4])
                    # remove unmapped reads
                    if len(bit_flag) < 7 or bit_flag[-3] == '1':
                        continue
                    # determine read is from read 1 or read 2
                    if len(bit_flag) >= 7 and bit_flag[-7] == '1':
                        sam_seq = str(index_r1[read_id].seq)
                    elif len(bit_flag) >= 7 and bit_flag[-7] == '0':
                        sam_seq = str(index_r2[read_id].seq)
                    # determine read is aligned to forward or reverse complement strand
                    if len(bit_flag) >= 7 and bit_flag[-5] == '1' and bit_flag[-7] == '1' and mapq >= 20:
                        sam_seq = str(Seq(sam_seq).reverse_complement())
                        # trim reads by cigar hard clipping 
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_rev.write(line.replace(line.split('\t')[9], sam_seq))
                    elif len(bit_flag) >= 8 and bit_flag[-5] == '1' and bit_flag[-8] == '1' and bit_flag[-2] == '1' and mapq >= 20:
                        sam_seq = str(Seq(sam_seq).reverse_complement())
                        # trim reads by cigar hard clipping 
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_fwd.write(line.replace(line.split('\t')[9], sam_seq))
                    elif len(bit_flag) >= 7 and bit_flag[-6] == '1' and bit_flag[-7] == '1' and bit_flag[-2] == '1' and mapq >= 20:
                        # trim reads by cigar hard clipping 
                        if cigar[-1] == 'H':
                            hc_right = int(re.split("([A-Z])",cigar)[-3])
                            sam_seq = sam_seq[0 : -hc_right]
                        if 'H' in cigar[0:-1]:
                            hc_left = int(cigar.split('H')[0])
                            sam_seq = sam_seq[hc_left:]
                        fout_fwd.write(line.replace(line.split('\t')[9], sam_seq))
                    elif len(bit_flag) >= 8 and bit_flag[-6] == '1' and bit_flag[-8] == '1' and bit_flag[-2] == '1' and mapq >= 20:
                        # trim reads by cigar hard clipping 
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