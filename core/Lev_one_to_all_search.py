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
import sys
import pymol
import pymol.cmd as cmd
from pymol import stored
import warnings
warnings.filterwarnings("ignore")



path= '/group/bioinf_protstr/Ballal/cath_trimmed_PDBs2/'


def main_pdb_seq(query_sequence):
        
        seq1_SS= query_sequence
        
        df_p2= df_seq 
        
        df_p2['SS_score']=-9
        df_p2['query_P']= sys.argv[1]

        for j in df_p2.index:
            seq1_len = len(seq1_SS)
            seq2_len = len(df_p2.loc[j,'SS_seq'])
            max_len = max(seq1_len, seq2_len)
            ss_score = abs(1-(lev(seq1_SS, df_p2.loc[j,'SS_seq'])/max_len)) #### SS_score calculation formula
            df_p2.loc[j,'SS_distance'] =   lev(seq1_SS, df_p2.loc[j,'SS_seq'])
            df_p2.loc[j,'SS_score'] =  ss_score
        df_p2 = df_p2[['id','SS_score','Class', 'Arch', 'Topol', 'Homol']]
        if int(sys.argv[2]) > 100:
            result = df_p2.nlargest(100, 'SS_score')
        else:    
            result = df_p2.nlargest(int(sys.argv[2]), 'SS_score')
        print('top results:\n',result.to_string(index=False)) #### show only top  SS-score values



if len(sys.argv) == 3:

    ##print('cmd entry:', sys.argv[1],type(sys.argv[1]))
    
    #file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_random_1000_speed_test_ForkPoolWorker-12.csv'
    file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv'
    df_seq= pd.read_csv(file_seq,skipinitialspace=True)
    df_seq=df_seq[~df_seq.SS_seq.isna()]
    df_seq=df_seq[~df_seq.AA_seq.isna()]
    df_seq.reset_index(inplace=True)
    df_seq =df_seq[['id','SS_seq']]
    df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
    #print(df_family.groupby('Homol')["Homol"].count().head(10))
    query_info = pd.DataFrame(df_family[df_family.id == sys.argv[1]])
    query_info = query_info[['id','Class', 'Arch', 'Topol', 'Homol']]
    print('query domain info:\n',query_info.to_string(index=False))
    df_seq =df_seq.merge(df_family, left_on='id', right_on='id', how='inner')

    p_path=path+sys.argv[1]+'.cif'
    pdb_id= sys.argv[1]
    if os.path.exists(p_path):
        cmd.load(p_path, pdb_id)
        sec_struc = []
        selection = "n. CA"
        cmd.iterate(selection, "sec_struc.append(ss)", space=locals())
        ss_seq = ''.join(['L' if a =='' else a for a in sec_struc])
        ##print('ss seq',ss_seq)
        cmd.delete(pdb_id)
        start = time.time()
        main_pdb_seq(ss_seq)
        print("lev time",time.time() - start)

else:
   
   print('two parameter is excepted') 
    

