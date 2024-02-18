
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

path_pdbs = '/group/bioinf_protstr/All_alphaFoldSS_Jun2023/codes/v2/v3/v4/v5/results/selected_represent_trimmed/'

no_threads = 1

path_tostore='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/'


df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/rad52_motif_list.csv')
df.reset_index(inplace=True)
# # initialize list elements
# data = ['A0A7S0JX89']
  
# # Create the pandas DataFrame with column name is provided explicitly
# df = pd.DataFrame(data, columns=['id'])


def main(df_merged_child):
    current_process=multiprocessing.current_process().name
    for i in df_merged_child.index:
        #p_path=path+df_merged_child.loc[i,'id']+'.cif'
        p_path=path_pdbs+df_merged_child.loc[i,'id']+'.pdb'
        pdb_id= df_merged_child.loc[i,'id']
        if os.path.exists(p_path):
            cmd.load(p_path, pdb_id)
            seq= cmd.get_fastastr(pdb_id)
            df_merged_child.loc[i,'AA_seq']= ''.join([a for a in seq.split('\n')[1:] if a!=''] )
            sec_struc = []
            selection = "n. CA"
            cmd.iterate(selection, "sec_struc.append(ss)", space=locals())
            df_merged_child.loc[i,'SS_seq'] = ''.join(['L' if a =='' else a for a in sec_struc])
            cmd.delete(pdb_id)
            sys.stdout.flush()
    # fileSave='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_'+current_process+'.csv'
    fileSave= path_tostore + current_process+'.csv'
    print(fileSave)
    df_merged_child.to_csv(fileSave, index=False)




#pool = multiprocessing.Pool(no_threads)



groups = np.array_split(df, no_threads)  # int(no_batches/no_threads)
pool = multiprocessing.Pool(no_threads)
#start = time.time()
pool.map(main, groups)
#print("tm time",time.time() - start)





exit()
#### file name in folder to csv

import os
#file_list = os.listdir('/group/bioinf_protstr/All_alphaFoldSS_Jun2023/codes/v2/v3/v4/v5/results/selected_represent_trimmed/')


file_list = ['.'.join(x.split('.')[:-1]) for x in os.listdir("/group/bioinf_protstr/All_alphaFoldSS_Jun2023/codes/v2/v3/v4/v5/results/selected_represent_trimmed/") 
             if os.path.isfile(os.path.join('/group/bioinf_protstr/All_alphaFoldSS_Jun2023/codes/v2/v3/v4/v5/results/selected_represent_trimmed/', x))]    
df = pd.DataFrame(file_list, columns=["id"])
df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/rad52_motif_list.csv')



### 
with open("/group/bioinf_protstr/To_ballal/orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/non-redundant-data-sets/cath-dataset-nonredundant-S40.list") as file:
    lines = [line.strip() for line in file]

df = pd.DataFrame(lines, columns=["id"])
df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_non_redundant_pdbs_s40_list.csv')
