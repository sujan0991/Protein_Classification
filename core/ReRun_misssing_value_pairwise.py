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

#######
#df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/ssap_ss_tm_local_results.csv')

# df = df[df['custom_local_idnt_q_s'].isna()]

# df = df[['Protein_q', 'Protein_s']]

# df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/ssap_local_missing_values.csv', index= False)


#######


#### this script is not complete ................

no_threads = 30
b = 0
no_batches=1

scores = [2, -1, -3, -2]

df_seq = pd.read_csv('/group/bioinf_protstr/toAli/stella_SS_seq_reduced.csv')
df_seq.reset_index(inplace=True)
df_seq =df_seq[['domain_id','sequence']]
path_tostore='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_DSSP_scope/'



missing_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/ssap_local_missing_values.csv')

missing_df = missing_df[:100]


#p_combinations_CATH_s20_LEV
def main(group_indices):
    c=scores
    for i in group_indices:
        df_p1= missing_df.loc[i]
        p1_name=df_p1.Protein_q
        seq_df1 = df_seq[df_seq["domain_id"] == p1_name]
        seq1_SS= seq_df1.sequence
        p2_name=df_p1.Protein_s
        seq_df2 = df_seq[df_seq["domain_id"] == p2_name]
        seq2_SS= seq_df2.sequence
        seq_df2['custom_local_idnt_p1_p2']=-9
        seq_df2['custom_local_idnt_p2_p1']=-9
        seq_df2['query_P']= p1_name
        seq_df2['query_len']=len(seq1_SS)
        seq_df2['subject_len']= len(seq2_SS)
        custom_local_alignments_p1_p2 = pairwise2.align.localms(seq1_SS, seq2_SS, c[0], c[1], c[2], c[3])
        if len(custom_local_alignments_p1_p2) > 0:
                best_custom_local_alignment_p1_p2 = max(custom_local_alignments_p1_p2, key=lambda alignment: alignment[2])
        else:
                best_custom_local_alignment_p1_p2 = [0, 0, 0]
        best_custom_local_alignment_ident_p1_p2 = best_custom_local_alignment_p1_p2[2] / max(len(seq1_SS), len(seq2_SS))
        #seq_df2['custom_local_idnt_p1_p2'] = best_custom_local_alignment_ident_p1_p2
        custom_local_alignments_p2_p1 = pairwise2.align.localms(seq2_SS, seq1_SS, c[0], c[1], c[2], c[3])
        if len(custom_local_alignments_p2_p1) > 0:
                best_custom_local_alignment_p2_p1 = max(custom_local_alignments_p2_p1, key=lambda alignment: alignment[2])
        else:
                best_custom_local_alignment_p2_p1 = [0, 0, 0]
        best_custom_local_alignment_ident_p2_p1 = best_custom_local_alignment_p2_p1[2] / max(len(seq1_SS), len(seq2_SS))
        #seq_df2['custom_local_idnt_p2_p1'] = best_custom_local_alignment_ident_p2_p1   
        print('p1_name,p2_name.....',i,p1_name,p2_name,best_custom_local_alignment_ident_p1_p2)
        sys.stdout.flush()
        # seq_df2 = seq_df2[['domain_id', 'query_P','query_len','subject_len','custom_local_idnt_p1_p2','custom_local_idnt_p2_p1']]
        # seq_df2.to_csv(path_tostore+'queryP_'+p1_name + 'batch_'+ str(b)+'.csv', index=False)
        
      


pool = multiprocessing.Pool(no_threads)   
P1= missing_df.index.to_list()
groups_indicesGlobal = np.array_split(P1, no_batches)  # int(no_batches/no_threads)
groups_indices = np.array_split(groups_indicesGlobal[b], no_threads)  # int(no_batches/no_threads)
pool.map(main, groups_indices)




