import numpy as np
import pandas as pd 
import glob
import re
import os
import pymol
import pymol.cmd as cmd



# ss_df = pd.read_csv("all_domain_sequence_original_SS.csv")
# print("ss len",len(ss_df))

aa_df = pd.read_csv("all_domain_sequence_original_AA.csv")
print("aa len",len(aa_df))

path = r'/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project'
filelist = [file for file in os.listdir(path) if file.startswith('tm_path_all_domain_with_BM_batch')]

print("filelist len",len(filelist))

for i,file in enumerate(filelist):
    
    single_csv = pd.read_csv(file)
    domain_id = single_csv["id"]
    domain_id_list =  list(domain_id)
    
    print("domain_id_list len",len(domain_id_list))
    
    single_batch =  aa_df[aa_df.domain_id.isin(domain_id_list)]

    print("single_batch len",len(single_batch))

    single_batch.to_csv('sequence_original_AA_batch-{}.csv'.format(i),index=False)


