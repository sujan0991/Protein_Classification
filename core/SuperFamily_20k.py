

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
import glob
import itertools
import glob
import os
import matplotlib.pyplot as plt

# secondFilter=True

# path= '/group/bioinf_protstr/Ballal/cath_pdb_first_modeled_list/' # use your path
# all_files = glob.glob((path+'*.csv'))
# li = []
# for filename in all_files:
#     df = pd.read_csv(filename, index_col=None, header=0)
#     li.append(df)

# df_merged = pd.concat(li, axis=0, ignore_index=True)

# print(len(df_merged))
# df_merged= df_merged[df_merged['first_modeled']==1]
# df_merged=df_merged[~df_merged['selection'].str.contains(",")]

# df_merged['length']= df_merged['selection'].str.split(' resi ',expand=True)[1]
# df_merged['length']= abs((df_merged['length'].str.split('-',expand=True)[0]).astype(int) - (df_merged['length'].str.split('-',expand=True)[1]).astype(int))



# df_merged= df_merged[(df_merged['length']>50) & (df_merged['length']<250)]
# df_merged.length.plot(kind='kde')
# plt.xlim(min(df_merged.length), max(df_merged.length))
# plt.savefig(path+'length.png', format="png")
# plt.close()
# print(len(df_merged))
# df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned-2.csv')

# exit()

# path= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/' # use your path

# # merge with super family
# # samples with having at least two from the same super family
# df_family= pd.read_csv(path+'cath-domain-description-v4_3_0.csv')
# df_merged =df_merged.merge(df_family, left_on='id', right_on='id', how='inner')

#  #df2 = pd.read_csv("")

# def sample(chunk, rate=0.3):
#     n = max(int(len(chunk)*rate), 1)
#     return chunk.sample(n=n, replace=True, random_state=1)

# df_merged=df_merged.groupby('Homol', group_keys=False).apply(sample)

# len(df_merged.Homol.unique())


# if secondFilter==True:
#     mask = (df_merged.groupby('Homol')['id']
#             .transform(lambda x: x.mask(x==0)    # mask out the 0 values
#                                     .nunique()     # count the nunique
#                         )
#             .gt(1)                                 # 1 pair
#         )
#     df_merged=df_merged[mask]

#     len(df_merged.Homol.unique())
# print(len(df_merged.Homol.unique()), len(df_merged))

# df_merged =df_merged[['id','pdb_id','selection','length','Class','Arch','Topol','Homol']]
# path2= '/group/bioinf_protstr/Ballal/cath_trimmed_PDBs/'
# df_merged['pdb_id']=path2+df_merged['id']+'.pdb'
# df_merged=df_merged[['id','pdb_id']]

# df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned_sample_20K_secondFilter:'+str(secondFilter)+'.csv', index=False)





####

# df_merged =pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned_sample_20K_secondFilter:'+str(secondFilter)+'.csv')

# domain_id = df_merged["id"]
# p_list =  list(domain_id)

# print("p_list length",len(p_list))

# p_combinations  = list(itertools.combinations(p_list,2))

# print("combination length",len(p_combinations))
# #path2= '/group/bioinf_protstr/Ballal/cath_trimmed_PDBs/'
# path2= '/group/bioinf_protstr/Ballal/pdb_trim_all/'
# p_combinations=pd.DataFrame(p_combinations)
# p_combinations.columns=['id1', 'id2']
# p_combinations['pdb_id_1']=path2+p_combinations['id1']+'.pdb'
# p_combinations['pdb_id_2']=path2+p_combinations['id2']+'.pdb'
# p_combinations.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations35K_secondFilter:'+str(secondFilter)+'.csv', index=False)


###






######### rest of the data

# def split_df(df):
#      if len(df) % 2 != 0:  # Handling `df` with `odd` number of rows
#       df = df.iloc[:-1, :]
#      df1, df2 =  np.array_split(df, 2)
#      return df1, df2




# batch0_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned_sample_20K_secondFilter:False.csv')

# batch0_df =batch0_df[['id']]
# batch0_df = batch0_df['id'].tolist()

# all_safe_cath =pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned.csv')

# not_batch0_df = all_safe_cath[~all_safe_cath.id.isin(batch0_df)]


# # shuffle the DataFrame rows
# not_batch0_df = not_batch0_df.sample(frac = 1)

# batch1_df,batch2_df = split_df(not_batch0_df)

# path2= '/group/bioinf_protstr/Ballal/pdb_trim_all/'

# batch1_df['pdb_id']=path2+batch1_df['id']+'.pdb'
# batch1_df=batch1_df[['id','pdb_id']]


# batch1_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned_batch1.csv', index=False)





# batch2_df['pdb_id']=path2+batch2_df['id']+'.pdb'
# batch2_df=batch2_df[['id','pdb_id']]

# batch2_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned_batch2.csv', index=False)




####### before this the code was for the result using feress method



##### below this is for ck1,ck2,ck3

no_con_no_con_df=pd.read_csv('/group/bioinf_protstr/ali/ck1-ck2-ck3_no_conflict.csv')
no_con_no_data_df= pd.read_csv('/group/bioinf_protstr/ali/ck2-ck3_no_conflict_ck1_noData.csv')


df_merged = no_con_no_data_df.append(no_con_no_con_df, ignore_index=True)

df_cath=pd.read_csv("/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-selection-v4_3_0.csv")
df_cath.columns = ['index', 'id', 'PDB_Id', 'selection']

df_merged= pd.merge(df_merged, df_cath ,how='left', on="PDB_Id")
df_merged=df_merged[~df_merged['selection'].str.contains(",")]

df_merged['length']= df_merged['selection'].str.split(' resi ',expand=True)[1]
df_merged['length']= abs((df_merged['length'].str.split('-',expand=True)[0]).astype(int) - (df_merged['length'].str.split('-',expand=True)[1]).astype(int) + 1)



df_merged= df_merged[(df_merged['length']>50) & (df_merged['length']<250)]
df_merged = df_merged[['PDB_Id','id', 'selection', 'length']]
df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/safe_cath_list_after_ck1-2-3.csv')


##df[df.length < 50]


# /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned_sample_20K_secondFilter:False.csv




# df_cath=df_cath[~df_cath['selection'].str.contains(",")]

