import argparse
import datetime
import gzip as gz
import os
import tempfile
import time
import multiprocessing as mp
import numpy as np
import pandas as pd
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import multiprocessing
import sys


def get_dssp(dssp_dict):
    dssp = []
    i_aa, i_d = 0, 1
    for key in dssp_dict.keys():
        dssp.append(dssp_dict[key][i_d])
    return  np.array(dssp)



def main(df_merged_child):
    current_process=multiprocessing.current_process().name
    print('......',len(df_merged_child))
    result_df = pd.DataFrame(columns=['domain_id','dssp'])
    for i in df_merged_child.index:
        id= df_merged_child.loc[i,'id']
        ##print(pdb_files+id+'.pdb')
        if os.path.exists(pdb_files+id+'.PDB'):
            print('file exist')
            path = pdb_files+id+'.PDB'
            cols = ['id','dssp']
            try:
                dssp_dict, _ = dssp_dict_from_pdb_file(path, DSSP='mkdssp')
                # aa, dssp = get_dssp(dssp_dict,chain,selection)
                dssp = get_dssp(dssp_dict)
                print('got DSSP',len(dssp))
                df = pd.DataFrame([(id, "".join(dssp))], columns=cols)
                result_df = pd.concat([result_df, df], axis=0)
            except (RuntimeError, Exception):
                # raise Exception("DSSP failed to produce an output for "+ path.split('/')[-1])
                print(f"DSSP failed to produce an output for {id}")
    result_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/Rad52DSSP_seq.csv',index=False)



#pdb_files = '/group/bioinf_protstr/Ballal/scope_pdb/'
pdb_files = '/group/bioinf_protstr/Ballal/test_pdbs2/'
no_threads = 1
no_batches=1


data = ['A0A7S0JX89']
  
# Create the pandas DataFrame with column name is provided explicitly
df_merged = pd.DataFrame(data, columns=['id'])


pool = multiprocessing.Pool(no_threads)
#df_merged = pd.read_csv('/group/bioinf_protstr/article_briefings/data/scope/scope.csv')
# df_merged = df_merged[['index','scop']]
# df_merged.columns = ['id','scop']
df_merged.reset_index(inplace=True)
groups = np.array_split(df_merged, no_batches)  # int(no_batches/no_threads)
pool = multiprocessing.Pool(no_threads)

start = time.time()
pool.map(main, groups)
print("tm time",time.time() - start)



exit()

df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_seq.csv')
df['dssp'] = df['dssp'].str.replace('-','N')
df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_seq_with_N.csv',index=False)
