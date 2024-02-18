import pandas as pd
import sys
import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np
import multiprocessing




file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned.csv'
df_seq= pd.read_csv(file_seq,skipinitialspace=True)
df_seq.drop_duplicates(subset='id', keep="last", inplace=True)
# df_seq=df_seq[~df_seq.SS_seq.isna()]
# df_seq=df_seq[~df_seq.AA_seq.isna()]
# df_seq=df_seq[~df_seq.id.isna()]
df_seq.reset_index(inplace=True)


df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-domain-description-v4_3_0.csv')
print(len(df_family.Homol.unique()))

df_seq =df_seq.merge(df_family, left_on='id', right_on='id', how='inner')
df_seq = df_seq[['id','Homol']]
print(len(df_seq.Homol.unique()))

total_combination = (len(df_seq)*(len(df_seq) - 1))/2




n_by_family = df_seq.groupby("Homol")["Homol"].count()

superfamily_df = pd.Series.to_frame(n_by_family)
# superfamily_df.rename_axis('S.Family')

superfamily_df['S.Family'] = list(n_by_family.index)
superfamily_df.reset_index(drop=True, inplace=True)

superfamily_df.columns = ['count','S.Family']
superfamily_df= superfamily_df[['S.Family','count']]

superfamily_df['sameFam']=((superfamily_df['count'])*(superfamily_df['count']-1))/2

print(superfamily_df.sum()[2])  ###214,049,278

diff = total_combination - superfamily_df.sum()[2]

print(diff)

# #### create just combonations for all
# b=2 ##### 0......2
# no_threads = 100
# no_batches=3

# file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned.csv'
# df_seq= pd.read_csv(file_seq,skipinitialspace=True)
# df_seq.reset_index(inplace=True)
# df_seq =df_seq[['id']]

# P2= set(df_seq.index)

# path_tostore='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/All_combination/'

# def main(group_indices):
#     current_process=multiprocessing.current_process().name
#     counter = 0
#     for i in group_indices:
#         counter = counter + 1 
#         df_p1= df_seq.loc[i] # first protein
#         p_name=df_p1.id
#         P2_sub_index = [x for x in  P2 if x >i]
#         df_p2= df_seq[df_seq.index.isin(P2_sub_index)] # keep only the proteins of interest to compare with
#         df_p2['query_P']=p_name
#         print(current_process,counter,len(group_indices), len(df_p2))
#         for j in df_p2.index:
#             print(j)
#         df_p2 = df_p2[['id', 'query_P']]
#         print('df_p2......',df_p2)
#         df_p2.to_csv(path_tostore+'queryP_'+p_name +'batch_'+ str(b)+'.csv', index=False)
#         sys.stdout.flush()
        
        
# pool = multiprocessing.Pool(no_threads)   
# P1= df_seq.index.to_list()
# groups_indicesGlobal = np.array_split(P1, no_batches)  # int(no_batches/no_threads)
# groups_indices = np.array_split(groups_indicesGlobal[b], no_threads)  # int(no_batches/no_threads)
# pool.map(main, groups_indices)        

