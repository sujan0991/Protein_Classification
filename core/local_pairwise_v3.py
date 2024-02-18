import numpy as np
import matplotlib.pyplot as plt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os, time, sys
import re, glob, random, datetime
from collections import defaultdict
from statistics import median, mean, stdev
import csv
import sys
import pandas as pd
from multiprocessing import Pool
import psutil
from pathlib import Path
import multiprocessing
from Levenshtein import distance as lev


match_dic = {}
match_dic[('S', 'S')] = 1
match_dic[('H', 'H')] = 1
match_dic[('L', 'L')] = 1
match_dic[('S', 'H')] = -1
match_dic[('H', 'S')] = -1
match_dic[('S', 'L')] = -0.5
match_dic[('L', 'S')] = -0.5
match_dic[('H', 'L')] = -0.5
match_dic[('L', 'H')] = -0.5




# sequences_df = pd.read_csv('/group/bioinf_protstr/toAli/stella_SS_seq_reduced.csv', usecols=[1,2])
#sequences_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_seq_reduced_ForkPoolWorker-1.csv')
## id,AA_seq,SS_seq

sequences_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_s20_DSSP_sequence_with_N.csv')

## id,dssp


def align_sequences_parallel(args):
    output_folder_path= '/group/bioinf/Ballal/p_com_local_s20_dssp/'
    column_name='dssp'
    query_string_full, query_id,sequences_df = args
    index_domain = sequences_df.index[sequences_df.id==query_id][0] # the index of current protein
    sequences_df=sequences_df[sequences_df.index>index_domain] # to avoid duplicated results
    sequences_df.index = sequences_df['id']
    sequences_df['query_id'] = query_id
    ##
    sequences_df['local_sim']= -9
    sequences_df['length_query']= -9
    sequences_df['length_target']= -9
    for i in sequences_df.index:
        sequence_full = sequences_df.loc[i, column_name]
        alignment = pairwise2.align.localds(sequence_full, query_string_full, match_dic,-1, -0.5, one_alignment_only=True)
        sequences_df.loc[i,'local_sim']=alignment[0][2]
        sequences_df.loc[i,'local_alignLength']=alignment[0][4]-alignment[0][3]
        sequences_df.loc[i,'length_query']=len(query_string_full)
        sequences_df.loc[i,'length_target']=len(sequence_full) 
    # Save the aligned DataFrame to a CSV file in the output folder
    sequences_df[['query_id', 'id','local_sim', 'local_alignLength','length_query', 'length_target']].to_csv(os.path.join(output_folder_path, 'aligned_sequences_' + str(query_id) + '.csv'), index=False)



# Number of processes in the pool

# output_folder_path= '/group/bioinf_protstr/toAli/local_pairswise_ssaps/'
# os.makedirs(output_folder_path, exist_ok=True)

align_args = [(qr_str,qr_id,sequences_df.copy()) for qr_str, qr_id in zip(sequences_df.dssp, sequences_df.id)]

## for batches
#align_args = [(qr_str, qr_id,sequences_df.copy()) for qr_str, qr_id in zip(sequences_df.sequence[9000:len(sequences_df)], sequences_df.domain_id[9000:len(sequences_df)])]

# Create a pool of processes and apply align_sequences_parallel to each argument tuple
with Pool(processes=240) as pool:
        pool.map(align_sequences_parallel, align_args)



