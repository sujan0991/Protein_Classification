import pandas as pd
import sys
import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np
import multiprocessing
import sys
import glob
import itertools
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns

length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv', usecols=['id','AA_seq','SS_seq'])
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')
length_df = length_df[['id','ss_length']]

df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')

### we need this for analisis

print(len(length_df.Class.unique()))
c_list = length_df.Class.value_counts().tolist()
same_class = 0
for i in c_list:
    print(i)
    c = (i*(i-1))/2
    same_class += c
print(same_class)    

print(len(length_df.Arch.unique()))
a_list = length_df.Arch.value_counts().tolist()
same_Arch = 0
for i in a_list:
    print(i)
    c = (i*(i-1))/2
    same_Arch += c
print(same_Arch)

print(len(length_df.Topol.unique()))
t_list = length_df.Topol.value_counts().tolist()
same_Topol = 0
for i in t_list:
    print(i)
    c = (i*(i-1))/2
    same_Topol += c
print(same_Topol)

print(len(length_df.Homol.unique()))
h_list = length_df.Homol.value_counts().tolist()

same_homol = 0
for i in h_list:
    print(i)
    c = (i*(i-1))/2
    same_homol += c
print(same_homol)


