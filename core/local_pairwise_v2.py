
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

## for customize
def perform_alignment(query_string, target_string, c=[2, -1, -3, -2]):
    #longer_seq_length = max(len(query_string),len(target_string))
    shorter_seq_length = min(len(query_string),len(target_string))
    local_alignments = pairwise2.align.localms(query_string, target_string, c[0], c[1], c[2], c[3], one_alignment_only=True)
    query_aln_string = local_alignments[0][0]
    target_aln_string = local_alignments[0][1]
    if len(local_alignments) > 0:local_alignments = local_alignments[0][2]
    else: local_alignments =0
    local_ident = local_alignments / (shorter_seq_length*c[0]) ## c[0] = match reward
    return local_alignments, local_ident, query_aln_string, target_aln_string

## for default
def perform_alignment(query_string, target_string):
    #longer_seq_length = max(len(query_string),len(target_string))
    shorter_seq_length = min(len(query_string),len(target_string))
    local_alignments = pairwise2.align.localxx(query_string, target_string, one_alignment_only=True)
    query_aln_string = local_alignments[0][0]
    target_aln_string = local_alignments[0][1]
    if len(local_alignments) > 0:local_alignments = local_alignments[0][2]
    else: local_alignments =0
    local_ident = local_alignments / (shorter_seq_length)
    return local_alignments, local_ident, query_aln_string, target_aln_string

def perform_alignment_global(query_string, target_string):
    longer_seq_length = max(len(query_string),len(target_string))
    #shorter_seq_length = min(len(query_string),len(target_string))
    local_alignments = pairwise2.align.globalxx(query_string, target_string, one_alignment_only=True)
    query_aln_string = local_alignments[0][0]
    target_aln_string = local_alignments[0][1]
    if len(local_alignments) > 0:local_alignments = local_alignments[0][2]
    else: local_alignments =0
    local_ident = local_alignments / (longer_seq_length)
    # lev_score = lev(s1, s2)
    # ss_score = abs(1 - lev_score / longer_seq_length)
    # print(ss_score)
    return local_alignments, local_ident, query_aln_string, target_aln_string



# sequences_df = pd.read_csv('/group/bioinf_protstr/toAli/stella_SS_seq_reduced.csv', usecols=[1,2])
sequences_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_seq_reduced_ForkPoolWorker-1.csv')
## id,AA_seq,SS_seq


def align_sequences_parallel(args):
    output_folder_path= '/group/bioinf/Ballal/p_com_local_s20/'
    column_name='SS_seq'
    column_name_reduced='SS_seq_reduced'
    query_string_full,query_string_reduced, query_id,sequences_df = args
    index_domain = sequences_df.index[sequences_df.id==query_id][0] # the index of current protein
    sequences_df=sequences_df[sequences_df.index>index_domain] # to avoid duplicated results
    sequences_df.index = sequences_df['id']
    sequences_df['query_id'] = query_id
    ##
    sequences_df['global_alignments_full'] = -9
    sequences_df['global_ident_full'] = -9
    sequences_df['query_aln_string_full'] = -9
    sequences_df['target_aln_string_full'] = -9
    ##
    sequences_df['global_alignments_reduced'] = -9
    sequences_df['global_ident_reduced'] = -9
    sequences_df['query_aln_string_reduced'] = -9
    sequences_df['target_aln_string_reduced'] = -9
    ##
    for i in sequences_df.index:
        sequence_full = sequences_df.loc[i, column_name]
        local_alignments, local_ident,query_aln_string_full, target_aln_string_full = perform_alignment_global(query_string_full, sequence_full)
        sequences_df.loc[i, 'global_alignments_full'] = local_alignments
        sequences_df.loc[i, 'global_ident_full'] = round(local_ident, 3)
        sequences_df.loc[i, 'query_aln_string_full'] = query_aln_string_full
        sequences_df.loc[i, 'target_aln_string_full'] = target_aln_string_full
        ##reduced
        sequence_reduced = sequences_df.loc[i, column_name_reduced]
        local_alignments_reduced, local_ident_reduced,query_aln_string_reduced, target_aln_string_reduced = perform_alignment_global(query_string_reduced, sequence_reduced)
        sequences_df.loc[i, 'global_alignments_reduced'] = local_alignments_reduced
        sequences_df.loc[i, 'global_ident_reduced'] = round(local_ident_reduced, 3)
        sequences_df.loc[i, 'query_aln_string_reduced'] = query_aln_string_reduced
        sequences_df.loc[i, 'target_aln_string_reduced'] = target_aln_string_reduced
    # Save the aligned DataFrame to a CSV file in the output folder
    sequences_df[['query_id', 'id','global_alignments_full', 'global_ident_full','query_aln_string_full', 'target_aln_string_full', 'global_alignments_reduced', 'global_ident_reduced','query_aln_string_reduced', 'target_aln_string_reduced' ]].to_csv(os.path.join(output_folder_path, 'aligned_sequences_' + str(query_id) + '.csv'), index=False)



# Number of processes in the pool

align_args = [(qr_str,qr_str_reduced,qr_id,sequences_df.copy()) for qr_str,qr_str_reduced, qr_id in zip(sequences_df.SS_seq,sequences_df.SS_seq_reduced, sequences_df.id)]

## for batches
#align_args = [(qr_str, qr_id,sequences_df.copy()) for qr_str, qr_id in zip(sequences_df.sequence[9000:len(sequences_df)], sequences_df.domain_id[9000:len(sequences_df)])]

# Create a pool of processes and apply align_sequences_parallel to each argument tuple
with Pool(processes=240) as pool:
        pool.map(align_sequences_parallel, align_args)




