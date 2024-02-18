
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

#### not working, problem in pdb file

def main(df_merged_child):
    current_process=multiprocessing.current_process().name
    for i in df_merged_child.index:
        p_path=path+df_merged_child.loc[i,'id']+'.pdb'
        print('path',p_path)
        pdb_id= df_merged_child.loc[i,'id']
        if os.path.exists(p_path):
            cmd.load(p_path, pdb_id)
            seq= cmd.get_fastastr(pdb_id)
            df_merged_child.loc[i,'AA_seq']= ''.join([a for a in seq.split('\n')[1:] if a!=''] )
            sec_struc = []
            selection = "n. CA"
            cmd.iterate(selection, "sec_struc.append(ss)", space=locals())
            df_merged_child.loc[i,'SS_seq'] = ''.join(['L' if a =='' else a for a in sec_struc])
            print(sec_struc)
            cmd.delete(pdb_id)
            sys.stdout.flush()
    # fileSave='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_'+current_process+'.csv'
    fileSave='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_'+current_process+'.csv'
    print(fileSave)
    df_merged_child.to_csv(fileSave, index=False)



path= '/group/bioinf_protstr/Ballal/scope_pdb/'

no_threads = 1
df = pd.read_csv('/group/bioinf_protstr/article_briefings/data/scope/scope.csv')
df = df[['index','scop']]
df.columns = ['id','scop']

df.reset_index(inplace=True)
groups = np.array_split(df, no_threads)  # int(no_batches/no_threads)
pool = multiprocessing.Pool(no_threads)
pool.map(main, groups)
