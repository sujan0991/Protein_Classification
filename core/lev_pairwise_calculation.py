# Imports

from multiprocessing import Pool, Manager
import pandas as pd
import os,time
from statistics import median,mean,stdev    # mean and median AlphaFold confidence
import matplotlib.pyplot as plt             # plot histograms etc.
import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
import glob
import multiprocessing
from Levenshtein import distance as lev
import psutil


no_threads = 60
b = 2
no_batches=3

scores = [2, -1, -3, -2]

##/group/bioinf_protstr/toAli/cath_SS_seq_reduced.csv      id,SS_seq,SuperFamily,SS_seq_reduced
##/group/bioinf_protstr/toAli/stella_SS_seq_reduced.csv              domain_id,sequence,SS_seq_reduced
df_seq = pd.read_csv('/group/bioinf_protstr/toAli/stella_SS_seq_reduced.csv')
df_seq.reset_index(inplace=True)
df_seq =df_seq[['domain_id','sequence']]
P2= set(df_seq.index)
path_tostore='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_DSSP_scope/'

#p_combinations_CATH_s20_LEV
def main(group_indices):
    c=scores
    for i in group_indices:
        df_p1= df_seq.loc[i] # first protein
        p_name=df_p1.domain_id
        seq1_SS= df_p1.sequence
        P2_sub_index = [x for x in  P2 if x >i]
        df_p2= df_seq[df_seq.index.isin(P2_sub_index)] # keep only the proteins of interest to compare with
        df_p2['custom_local_idnt_p1_p2']=-9
        df_p2['custom_local_idnt_p2_p1']=-9
        df_p2['query_P']=p_name
        df_p2['query_len']=len(seq1_SS)
        df_p2['subject_len']=-9
        for j in df_p2.index:
            seq2_SS = df_p2.loc[j,'sequence']
            df_p2.loc[j,'subject_len'] = len(seq2_SS)
            ## custom_local_alignments_p1_p2
            custom_local_alignments_p1_p2 = pairwise2.align.localms(seq1_SS, seq2_SS, c[0], c[1], c[2], c[3])
            if len(custom_local_alignments_p1_p2) > 0:
                best_custom_local_alignment_p1_p2 = max(custom_local_alignments_p1_p2, key=lambda alignment: alignment[2])
            else:
                best_custom_local_alignment_p1_p2 = [0, 0, 0]
            best_custom_local_alignment_ident_p1_p2 = best_custom_local_alignment_p1_p2[2] / max(len(seq1_SS), len(seq2_SS))
            df_p2.loc[j,'custom_local_idnt_p1_p2'] = best_custom_local_alignment_ident_p1_p2
            ## custom_local_alignments_p2_p1
            custom_local_alignments_p2_p1 = pairwise2.align.localms(seq2_SS, seq1_SS, c[0], c[1], c[2], c[3])
            if len(custom_local_alignments_p2_p1) > 0:
                best_custom_local_alignment_p2_p1 = max(custom_local_alignments_p2_p1, key=lambda alignment: alignment[2])
            else:
                best_custom_local_alignment_p2_p1 = [0, 0, 0]
            best_custom_local_alignment_ident_p2_p1 = best_custom_local_alignment_p2_p1[2] / max(len(seq1_SS), len(seq2_SS))
            df_p2.loc[j,'custom_local_idnt_p2_p1'] = best_custom_local_alignment_ident_p2_p1    
        df_p2 = df_p2[['domain_id', 'query_P','query_len','subject_len','custom_local_idnt_p1_p2','custom_local_idnt_p2_p1']]
        df_p2.to_csv(path_tostore+'queryP_'+p_name + 'batch_'+ str(b)+'.csv', index=False)
        
      

pool = multiprocessing.Pool(no_threads)   
P1= df_seq.index.to_list()
groups_indicesGlobal = np.array_split(P1, no_batches)  # int(no_batches/no_threads)
groups_indices = np.array_split(groups_indicesGlobal[b], no_threads)  # int(no_batches/no_threads)
pool.map(main, groups_indices)

## local gives similarity as gap opening is -




