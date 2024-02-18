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


# H_alph = ['G', 'H', 'I']
# S_alph = ['E', 'B']
# L_alph = ['S', 'T', '-']


def get_dssp(dssp_dict,chain,selection):

    first_selection = selection.split('-')[0]
    last_selection = selection.split('-')[1]

    ##   ('A', (' ', 187, ' '))
    # ## ('A', (' ', 1, ' '))
    firsr_key = ('%s'%chain, (' ', int(first_selection), ' '))
    last_key = ('%s'%chain, (' ', int(last_selection), ' '))
    print(firsr_key)
    print(last_key)

    kys = list(dssp_dict.keys())
    print('.........................',kys)
    kys = kys[kys.index(firsr_key):kys.index(last_key)+1]
    
    newDict = {k:dssp_dict[k] for k in kys}
    
    dssp_dict = newDict

    
    # aa = []
    dssp = []
    i_aa, i_d = 0, 1

    for key in dssp_dict.keys():
        #print('dssp....',dssp_dict)
        ##aa.append(dssp_dict[key][i_aa])
        dssp.append(dssp_dict[key][i_d])
    # return np.array(aa), np.array(dssp)
    return  np.array(dssp)



    

# path = '/group/bioinf_protstr/Ballal/1a0f.cif'

# print('??????????????????????????',_row_fn('1a7x','A','1-107',path))

# exit()

def main(df_merged_child):
    current_process=multiprocessing.current_process().name
    print('......',len(df_merged_child))

    result_df = pd.DataFrame(columns=['domain_id','dssp'])

    for i in df_merged_child.index:
        pdb_id=df_merged_child.loc[i,'PDB_Id']
        id= df_merged_child.loc[i,'id']
        selection = df_merged_child.loc[i,'selection']

        chain = selection.split(' and ')[0]
        selection = selection.split(' and ')[1]
        chain = chain[-1]
        selection = selection.split(' ')[1]
        
        if os.path.exists(pdb_files+pdb_id+'.cif'):
            print('file exist',pdb_id)
            path = pdb_files+pdb_id+'.cif'
            
            cols = ['domain_id','dssp']

            # try:
            dssp_dict, _ = dssp_dict_from_pdb_file(path, DSSP='mkdssp')
                # aa, dssp = get_dssp(dssp_dict,chain,selection)
            dssp = get_dssp(dssp_dict,chain,selection)
            print('got DSSP',len(dssp))
            df = pd.DataFrame([(id, "".join(dssp))], columns=cols)
            result_df = pd.concat([result_df, df], axis=0)
            # except (RuntimeError, Exception):
            #     # raise Exception("DSSP failed to produce an output for "+ path.split('/')[-1])
            #     print(f"DSSP failed to produce an output for {id}")

                ##return pd.DataFrame([(id, "")], columns=cols)
        
    result_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_DSSP_seq_after_ck1-2-3.csv',index=False)



pdb_files = '/group/bioinf_protstr/Ballal/all_pdbs/'
no_threads = 1
no_batches=1

pool = multiprocessing.Pool(no_threads)
df_merged=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/safe_cath_list_after_ck1-2-3.csv')
df_merged.reset_index(inplace=True)
groups = np.array_split(df_merged, no_batches)  # int(no_batches/no_threads)
pool = multiprocessing.Pool(no_threads)
pool.map(main, groups)





#### replce '-' by N

df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_DSSP_seq_after_ck1-2-3.csv')
df['dssp'] = df['dssp'].str.replace('-','N')
df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_DSSP_seq_with_N_after_ck1-2-3.csv',index=False)






