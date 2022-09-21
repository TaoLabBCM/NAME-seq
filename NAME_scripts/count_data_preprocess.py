#!/usr/bin/python

import sys
import re
import pandas as pd
import numpy as np

native_readcount_fwd = pd.read_csv(sys.argv[1], sep="\t", header=None, skiprows = [0,1,2], usecols=range(5), index_col = 0, names=['Position','A','C','G','T'])
native_readcount_rev = pd.read_csv(sys.argv[2], sep="\t", header=None, skiprows = [0,1,2], usecols=range(5), index_col = 0, names=['Position','A','C','G','T'])
pcr_readcount_fwd = pd.read_csv(sys.argv[3], sep="\t", header=None, skiprows = [0,1,2], usecols=range(5), index_col = 0, names=['Position','A','C','G','T'])
pcr_readcount_rev = pd.read_csv(sys.argv[4], sep="\t", header=None, skiprows = [0,1,2], usecols=range(5), index_col = 0, names=['Position','A','C','G','T'])
ref_fasta = sys.argv[5]
output_file = sys.argv[6]
with open(ref_fasta, 'r') as ref:
    ref_seq=''
    for line in ref:
        if line[0] == '>':
            chrom_name = line.rstrip().split(' ')[0][1:]
        if line[0] != '>':
            ref_seq += line.rstrip()
nrow = len(ref_seq)
if len(native_readcount_fwd) < nrow:
    missing_row_index = list(set(range(1,nrow + 1)) - set(native_readcount_fwd.index.tolist()))
    df_with_missing_row = pd.DataFrame(np.zeros((len(missing_row_index),len(native_readcount_fwd.columns.tolist()))), 
                                       index = missing_row_index, 
                                       columns = native_readcount_fwd.columns.tolist())
    native_readcount_fwd = pd.concat([native_readcount_fwd, df_with_missing_row]).sort_index()
if len(native_readcount_rev) < nrow:
    missing_row_index = list(set(range(1,nrow + 1)) - set(native_readcount_rev.index.tolist()))
    df_with_missing_row = pd.DataFrame(np.zeros((len(missing_row_index),len(native_readcount_rev.columns.tolist()))), 
                                       index = missing_row_index, 
                                       columns = native_readcount_rev.columns.tolist())
    native_readcount_rev = pd.concat([native_readcount_rev, df_with_missing_row]).sort_index()
if len(pcr_readcount_fwd) < nrow:
    missing_row_index = list(set(range(1,nrow + 1)) - set(pcr_readcount_fwd.index.tolist()))
    df_with_missing_row = pd.DataFrame(np.zeros((len(missing_row_index),len(pcr_readcount_fwd.columns.tolist()))), 
                                       index = missing_row_index, 
                                       columns = pcr_readcount_fwd.columns.tolist())
    pcr_readcount_fwd = pd.concat([pcr_readcount_fwd, df_with_missing_row]).sort_index()
if len(pcr_readcount_rev) < nrow:
    missing_row_index = list(set(range(1,nrow + 1)) - set(pcr_readcount_rev.index.tolist()))
    df_with_missing_row = pd.DataFrame(np.zeros((len(missing_row_index),len(pcr_readcount_rev.columns.tolist()))), 
                                       index = missing_row_index, 
                                       columns = pcr_readcount_rev.columns.tolist())
    pcr_readcount_rev = pd.concat([pcr_readcount_rev, df_with_missing_row]).sort_index()
reverse_complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
# chrom = native_readcount_fwd.iat[0,0]
data = []
chrom = chrom_name
for i in range(nrow):
    # forward strand preprocess
    start_loc = i
    end_loc = start_loc + 1
    strand = '+'
    ref = ref_seq[i]
    # get original basa count
    native_A = native_readcount_fwd.iat[i,0]
    native_C = native_readcount_fwd.iat[i,1]
    native_G = native_readcount_fwd.iat[i,2]
    native_T = native_readcount_fwd.iat[i,3]
#         native_err_count = sum(sorted(set([native_A, native_C, native_G, native_T]))[0:2])/2
#         # calculate corrected base count
#         native_A = max([native_A - native_err_count, 0.0])
#         native_C = max([native_C - native_err_count, 0.0])
#         native_G = max([native_G - native_err_count, 0.0])
#         native_T = max([native_T - native_err_count, 0.0])
    native_cor = sum([native_A, native_C, native_G, native_T])
    pcr_A = pcr_readcount_fwd.iat[i,0]
    pcr_C = pcr_readcount_fwd.iat[i,1]
    pcr_G = pcr_readcount_fwd.iat[i,2]
    pcr_T = pcr_readcount_fwd.iat[i,3]
#         pcr_err_count = sum(sorted(set([pcr_A, pcr_C, pcr_G, pcr_T]))[0:2])/2
#         pcr_A = max([pcr_A - pcr_err_count, 0.0])
#         pcr_C = max([pcr_C - pcr_err_count, 0.0])
#         pcr_G = max([pcr_G - pcr_err_count, 0.0])
#         pcr_T = max([pcr_T - pcr_err_count, 0.0])
    pcr_cor = sum([pcr_A, pcr_C, pcr_G, pcr_T])
    data.append([chrom, start_loc, end_loc, strand, ref, native_A, native_C, native_G, native_T, pcr_A, pcr_C, pcr_G, pcr_T, native_cor, pcr_cor])
    # reverse strand preprocess
    strand = '-'
    ref =  reverse_complement[ref_seq[i]]
    # get original basa count
    native_T = native_readcount_rev.iat[i,0]
    native_G = native_readcount_rev.iat[i,1]
    native_C = native_readcount_rev.iat[i,2]
    native_A = native_readcount_rev.iat[i,3]
#         native_err_count = sum(sorted(set([native_A, native_C, native_G, native_T]))[0:2])/2
#         # calculate corrected base count
#         native_A = max([native_A - native_err_count, 0.0])
#         native_C = max([native_C - native_err_count, 0.0])
#         native_G = max([native_G - native_err_count, 0.0])
#         native_T = max([native_T - native_err_count, 0.0])
    native_cor = sum([native_A, native_C, native_G, native_T])
    pcr_T = pcr_readcount_rev.iat[i,0]
    pcr_G = pcr_readcount_rev.iat[i,1]
    pcr_C = pcr_readcount_rev.iat[i,2]
    pcr_A = pcr_readcount_rev.iat[i,3]
#         pcr_err_count = sum(sorted(set([pcr_A, pcr_C, pcr_G, pcr_T]))[0:2])/2
#         pcr_A = max([pcr_A - pcr_err_count, 0.0])
#         pcr_C = max([pcr_C - pcr_err_count, 0.0])
#         pcr_G = max([pcr_G - pcr_err_count, 0.0])
#         pcr_T = max([pcr_T - pcr_err_count, 0.0])
    pcr_cor = sum([pcr_A, pcr_C, pcr_G, pcr_T])
    data.append([chrom, start_loc, end_loc, strand, ref, native_A, native_C, native_G, native_T, pcr_A, pcr_C, pcr_G, pcr_T, native_cor, pcr_cor])
preprocessed_data = pd.DataFrame(data, columns = ['chrom', 'start_loc', 'end_loc', 'strand', 'ref',
                                         'native_A', 'native_C','native_G', 'native_T', 
                                         'pcr_A', 'pcr_C','pcr_G', 'pcr_T', 
                                         'native_cor', 'pcr_cor'])
preprocessed_data.to_csv(output_file, sep = '\t', index = False)
