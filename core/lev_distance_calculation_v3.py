
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
from Levenshtein import distance as lev
import warnings
warnings.filterwarnings("ignore")


b= 1
no_threads = 30
no_batches=2

path_tostore='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_CATH_s20_LEV/'


file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/cath_non_redundant_pdbs_s20_sequence_reduced_ForkPoolWorker-1.csv'
df_seq= pd.read_csv(file_seq)
df_seq.reset_index(inplace=True)
df_seq =df_seq[['id','SS_seq_reduced']]
P2= set(df_seq.index)



def main(group_indices):
    current_process=multiprocessing.current_process().name
    counter = 0
    for i in group_indices:
        counter = counter + 1 
        df_p1= df_seq.loc[i] # first protein
        p_name=df_p1.id
       # seq1_AA= df_p1.AA_seq
        seq1_SS= df_p1.SS_seq_reduced
        P2_sub_index = [x for x in  P2 if x >i]
        df_p2= df_seq[df_seq.index.isin(P2_sub_index)] # keep only the proteins of interest to compare with
        #df_p2['AA_distance']=-9
        df_p2['SS_distance']=-9
        df_p2['query_P']=p_name
        #df_p2['query_len']=len(seq1_SS)
        print(current_process,counter)
        for j in df_p2.index:
           # df_p2.loc[j,'AA_distance'] =   lev(seq1_AA, df_p2.loc[j,'AA_seq'])
            df_p2.loc[j,'SS_distance'] =   lev(seq1_SS, df_p2.loc[j,'SS_seq_reduced'])
            #df_p2['subject_len']= len(df_p2.loc[j,'SS_seq'])
        df_p2 = df_p2[['id', 'query_P','SS_distance']]
        df_p2.to_csv(path_tostore+'queryP_'+p_name +'batch_'+ str(b)+'.csv', index=False)
        sys.stdout.flush()



pool = multiprocessing.Pool(no_threads)   
P1= df_seq.index.to_list()
groups_indicesGlobal = np.array_split(P1, no_batches)  # int(no_batches/no_threads)
groups_indices = np.array_split(groups_indicesGlobal[b], no_threads)  # int(no_batches/no_threads)
pool.map(main, groups_indices)






