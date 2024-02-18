# lev distance calculation
import pandas as pd
import sys
import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np
import multiprocessing
import sys
import pymol
import pymol.cmd as cmd
from pymol import stored
from Levenshtein import distance as lev

b = 0

file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_ForkPoolWorker-2.csv'
df_seq= pd.read_csv(file_seq,skipinitialspace=True)
df_seq.drop_duplicates(subset='id', keep="last", inplace=True)
df_seq=df_seq[~df_seq.SS_seq.isna()]
df_seq=df_seq[~df_seq.AA_seq.isna()]
df_seq=df_seq[~df_seq.id.isna()]

import warnings
warnings.filterwarnings("ignore")

def main(p_combinations_child):
    current_process=multiprocessing.current_process().name
    counter = 0
    p_combinations_child['AA_distance']=-9
    p_combinations_child['SS_distance']=-9
    for i in p_combinations_child.index:
        counter = counter + 1 
        if (counter % 1000) == 0:print(current_process,counter,len(p_combinations_child))
        sys.stdout.flush()
        pdb_id1= p_combinations_child.loc[i,'id1']
        pdb_id2= p_combinations_child.loc[i,'id2']
        row = df_seq[df_seq.id==pdb_id1]
        AA_seq1 = row.AA_seq.item()
        SS_seq1 = row.SS_seq.item()
        row = df_seq[df_seq.id==pdb_id2]
        if not row.empty:
         AA_seq2 = row.AA_seq.item()
         SS_seq2 = row.SS_seq.item()
         p_combinations_child.loc[i,'AA_distance'] = lev(AA_seq1, AA_seq2)
         p_combinations_child.loc[i,'SS_distance'] = lev(SS_seq1, SS_seq2)

        else:
                        
         print('.........current_process',i,current_process,pdb_id1,pdb_id2)

   
    p_combinations_child.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_LEV_childrenFilter_False/'+current_process +'_'+ str(b)+'.csv', index=False)


no_threads = 30
no_batches=5
pool = multiprocessing.Pool(no_threads)
p_combinations= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations35K_secondFilter:False.csv', usecols=['id1','id2'])
p_combinations.reset_index(inplace=True)
groupsGlobal = np.array_split(p_combinations, no_batches)  # int(no_batches/no_threads)
p_combinations=0
groups = np.array_split(groupsGlobal[b], no_threads)  # int(no_batches/no_threads)
groupsGlobal=0
pool.map(main, groups)























