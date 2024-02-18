

# import psico.editing

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

def first_modeled(id):
    # find index of first modeled residue
    var_space = {'i0': []}
    cmd.iterate('first %s and present' % id, 'i0.append(resi)', space=var_space)
    try: i0 = int(var_space['i0'][0])
    except: i0=-99
    return i0


def main(df_cath_child):
    current_process=multiprocessing.current_process().name
    for i in df_cath_child.index:
        pdb_id=df_cath_child.loc[i,'pdb_id']
        chain=df_cath_child.loc[i,'chain']
        if os.path.exists(path+pdb_id+'.cif.gz'):
            cmd.load(path+pdb_id+'.cif.gz',pdb_id)
            selection = pdb_id+ ' and ' +chain
            selection_name = pdb_id+'-'+'chain'
            cmd.select(selection_name, selection)
            fmi= first_modeled(selection_name)
            df_cath_child.loc[i,'first_modeled']=int(fmi)
            cmd.delete(pdb_id)
            if (i % 1000) == 0:print(current_process,i,len(df_cath_child))
            sys.stdout.flush()
    df_cath_child.to_csv('/group/bioinf_protstr/Ballal/cath_pdb_first_modeled_list/'+current_process+'.csv')

   




no_threads = 30
no_batches=30
pool = multiprocessing.Pool(no_threads)
df_cath=pd.read_csv("/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-selection-v4_3_0.csv")
df_cath['chain']= df_cath['selection'].str.split(' and',expand=True)[0]
df_cath['first_modeled']=-8
path= '/group/bioinf_protstr/Ballal/ali/'

groups = np.array_split(df_cath, no_batches)  # int(no_batches/no_threads)
pool = multiprocessing.Pool(no_threads)
pool.map(main, groups)
