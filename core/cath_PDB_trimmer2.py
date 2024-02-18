
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



def main(df_merged_child):
    current_process=multiprocessing.current_process().name
    for i in df_merged_child.index:
        pdb_id=df_merged_child.loc[i,'PDB_Id']
        id= df_merged_child.loc[i,'id']
        selection=df_merged_child.loc[i,'selection']
        if os.path.exists(path+pdb_id+'.cif'):
            cmd.load(path+pdb_id+'.cif',pdb_id)
            selection = pdb_id+ ' and ' +selection
            cmd.create(id, selection)
            cmd.save(path2+id+'.cif', id)
            cmd.delete(pdb_id)
            cmd.delete(id)
            if (i % 100) == 0:print(current_process,i,len(df_merged_child))
            sys.stdout.flush()





path= '/group/bioinf_protstr/Ballal/all_pdbs/'
path2= '/group/bioinf_protstr/Ballal/cath_trimmed_PDBs2/'

no_threads = 20
no_batches=20
pool = multiprocessing.Pool(no_threads)
df_merged=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/safe_cath_list_after_ck1-2-3.csv')
df_merged.reset_index(inplace=True)
groups = np.array_split(df_merged, no_batches)  # int(no_batches/no_threads)
pool = multiprocessing.Pool(no_threads)
pool.map(main, groups)














