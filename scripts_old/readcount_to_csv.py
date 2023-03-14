#!/usr/bin/python

import sys
import re
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import csv

# get igvtools readcount output, reference genome, and csv output file
readcount_fwd = sys.argv[1]
readcount_rev = sys.argv[2]
ref_genome = sys.argv[3]
csv_data_output = sys.argv[4]

def wig_to_csv(wig_file_fwd, wig_file_rev, chrom, cor_threshold = 0):
    with open(wig_file_fwd,'r') as wig_fwd:
        with open(wig_file_rev,'r') as wig_rev:
            with open(csv_data_output, 'a') as csv_f:
                chrom_found = False
                line_output = []

                for line in wig_fwd:
                    if line.startswith('variable'):
                        if line.split(' ')[1].split('=')[1] == chrom:
                            chrom_found = True
                            continue
                    if chrom_found:
                        if line.startswith('variable'):
                            break
                        line_list = line.rstrip().split('\t')
                        ref = ref_fasta_dict[chrom][int(line_list[0])-1]
                        if ref == 'A':
                            cor = float(line_list[1]) + float(line_list[2]) + float(line_list[3]) + float(line_list[4])
                            if cor <= cor_threshold:
                                continue
                            AtoG = float(line_list[3])/cor
                            AtoT = float(line_list[4])/cor
                            line_output.append([chrom, line_list[0], ref, '+', line_list[1], line_list[2], line_list[3], line_list[4], str(cor), str(AtoG), str(AtoT)])
                write = csv.writer(csv_f)
                write.writerows(line_output)
                
                chrom_found = False
                line_output = []
                for line in wig_rev:
                    if line.startswith('variable'):
                        if line.split(' ')[1].split('=')[1] == chrom:
                            chrom_found = True
                            continue
                    if chrom_found:
                        if line.startswith('variable'):
                            break
                        line_list = line.rstrip().split('\t')
                        ref = ref_fasta_dict[chrom][int(line_list[0])-1]
                        if ref == 'T':
                            cor = float(line_list[1]) + float(line_list[2]) + float(line_list[3]) + float(line_list[4])
                            if cor <= cor_threshold:
                                continue
                            AtoG = float(line_list[2])/cor
                            AtoT = float(line_list[1])/cor
                            line_output.append([chrom, line_list[0], ref, '-', line_list[1], line_list[2], line_list[3], line_list[4], str(cor), str(AtoG), str(AtoT)])
                write = csv.writer(csv_f)
                write.writerows(line_output)
    return()

# Load reference genome 
ref_fasta_dict = {rec.id : rec.seq for rec in SeqIO.parse(ref_genome, "fasta")}
# Run the readcount_to_csv function
for chrom in ref_fasta_dict.keys():
    wig_to_csv(readcount_fwd, readcount_rev, chrom)