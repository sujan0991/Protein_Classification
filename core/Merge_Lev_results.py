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


## for ferras analisis

#path= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_all_LEV_childrenFilter_False_v3/' # use your path
path= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_all_LEV_after_ck1-2-3/'
all_files = glob.glob((path+'*.csv'))
print('all_files len',len(all_files))
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    #length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/merged_Cath_list_cleaned.csv')
    length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv')
    length_df = length_df[['id','length']]
    df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
    length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')
    length_df =length_df[['id','length','Homol']]
    df = pd.merge(df, length_df ,how='left', on="id")
    df = df[['id', 'query_P', 'AA_distance', 'SS_distance', 'length', 'Homol']]
    df.columns = ['domain1','id','AA_distance', 'SS_distance', 'domain1_length','Homol']
    df = pd.merge(df, length_df ,how='left', on="id")
    df.columns = ['domain1','domain2','AA_distance', 'SS_distance', 'domain1_length', 'Homol1','domain2_length', 'Homol2']
    df['max_seq_len'] = df[["domain1_length", "domain2_length"]].max(axis=1)
    df['lev%_SS'] = df['SS_distance']/ df['max_seq_len']
    #### x
    df['SS_score'] = abs(1-df["lev%_SS"]).astype(np.float16)
    df['lev%_AA'] = df['AA_distance']/ df['max_seq_len']
    #### x
    df['AA_score'] = abs(1-df["lev%_AA"]).astype(np.float16)
    df['cath_superFamily'] = np.where(df["Homol1"] == df["Homol2"], 1, 0)
    #df = df[['SS_score','cath_superFamily']]
    df = df[['domain1','domain2', 'cath_superFamily','SS_score','AA_score']]
    li.append(df)

df_merged = pd.concat(li, axis=0, ignore_index=True)

df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/lev_results_all_safe_cath_after_ck1-2-3_v2.csv')

print('len(df_merged)',len(df_merged))

exit()



####### check why some score are over 1

df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_all_LEV_after_ck1-2-3_v2/queryP_1d6rI00batch_0.csv')
length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv')
length_df = length_df[['id','length']]
df = pd.merge(df, length_df ,how='left', on="id")
df.columns = ['domain1','id','AA_distance', 'SS_distance', 'domain1_length']
df = pd.merge(df, length_df ,how='left', on="id")
df.columns = ['domain1','domain2','AA_distance', 'SS_distance', 'domain1_length','domain2_length']
df['max_seq_len'] = df[["domain1_length", "domain2_length"]].max(axis=1)
df['lev%_SS'] = df['SS_distance']/ df['max_seq_len']
df['SS_score'] = abs(1-df["lev%_SS"]).astype(np.float16)
high_ss = df.query('SS_score > 1')
test = df.loc[df['domain1'] == '5a96A00']
test2 = df_merged.loc[(df_merged['domain1'] == '5a96A00') & (df_merged['domain2'] == '1d6rI00')]
test3 = df.loc[df['domain1'] == '4tz3A01']
test4= df_merged.loc[(df_merged['domain1'] == '4tz3A01') & (df_merged['domain2'] == '1d6rI00')]
##########################

df_merged= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/lev_results_all_safe_cath_after_ck1-2-3_v2.csv')
#df_merged['SS_score'].max()
#low_ss_aa = df_merged.query('SS_score <= 1' and 'AA_score <= 1')
## df_merged = df_merged.query('SS_score > 1' and 'AA_score > 1')

same_homo = df_merged.query('cath_superFamily == 1') ['SS_score']
diff_homo = df_merged.query('cath_superFamily == 0') ['SS_score']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 0.5', above_point_5 ,'below point 0.5', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 0.5', above_point_5 ,'below point 0.5', below_point_5)
 
same_homo_df = pd.DataFrame({'SS_score':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'SS_score':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["SS_score"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_score_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_score_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()

#################




path= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_all_LEV_childrenFilter_False_v3/' # use your path
all_files = glob.glob((path+'*.csv'))
print('all_files len',len(all_files))
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    
    li.append(df)

df_merged = pd.concat(li, axis=0, ignore_index=True)

df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_results_all_safe_cath_childrenFilter_False_v3.csv')

print('len(df_merged)',len(df_merged))

exit()

df_merged = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_results_all_childrenFilter_False_v3.csv')### test set all
length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned.csv')
length_df = length_df[['id','length']]
df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')
length_df =length_df[['id','length','Homol']]
df_merged = pd.merge(df_merged, length_df ,how='left', on="id")
df_merged = df_merged[['id', 'query_P', 'AA_distance', 'SS_distance', 'length', 'Homol']]
df_merged.columns = ['domain1','id','AA_distance', 'SS_distance', 'domain1_length','Homol']
df_merged = pd.merge(df_merged, length_df ,how='left', on="id")
#print(df_merged.isnull().any())
#['domain1', 'id', 'AA_distance', 'SS_distance', 'domain1_length', 'Homol_x', 'length', 'Homol_y']
df_merged.columns = ['domain1','domain2','AA_distance', 'SS_distance', 'domain1_length', 'Homol1','domain2_length', 'Homol2']


##df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_results_all_with_superfamily_childrenFilter_False_v3.csv',index=False)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from scipy.stats import pearsonr



df_merged['max_seq_len'] = df_merged[["domain1_length", "domain2_length"]].max(axis=1)
df_merged['lev%_SS'] = df_merged['SS_distance']/ df_merged['max_seq_len']
#### x
df_merged['seq_idnt'] = abs(1-df_merged["lev%_SS"])

same_homo = df_merged.query('Homol1 == Homol2') ['seq_idnt']
diff_homo = df_merged.query('Homol1 != Homol2') ['seq_idnt']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 0.5', above_point_5 ,'below point 0.5', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 0.5', above_point_5 ,'below point 0.5', below_point_5)
 
# same_homo_df = pd.DataFrame({'seq_idnt_SS':same_homo.values, 'Catagory':'Same Superfamily'})
# diff_homo_df = pd.DataFrame({'seq_idnt_SS':diff_homo.values, 'Catagory':'Different Superfamily'})
 
# same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

# g2 = sns.violinplot( y=same_homo_not_same_homo_df["seq_idnt_SS"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# #g2.set(ylabel=None)
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/seq_idnt_SS_super_safe_cath2.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/seq_idnt_SS_super_safe_cath2.svg', format="svg")
# plt.close()
# same_homo_not_same_homo_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/same_different_superfamily.csv',index=False)

# df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_results_all_with_superfamily_childrenFilter_False_v3_onlyHomo.csv',index=False)
#####



##### AA


df_merged['lev%_AA'] = df_merged['AA_distance']/ df_merged['max_seq_len']
#### x
df_merged['AA_seq_idnt'] = abs(1-df_merged["lev%_AA"])

same_homo = df_merged.query('Homol1 == Homol2') ['AA_seq_idnt']
diff_homo = df_merged.query('Homol1 != Homol2') ['AA_seq_idnt']

below_point_3 = 0
above_point_3 = 0
for i, v in same_homo.items():
    if v > 0.3:
        above_point_3 = above_point_3 + 1
    else:
        below_point_3 = below_point_3 + 1   

print('lev AA for same superfamily:total',len(same_homo), 'above point 0.3', above_point_3 ,'below point 0.3', below_point_3)


below_point_3 = 0
above_point_3 = 0
for i, v in diff_homo.items():
    if v > 0.3:
        above_point_3 = above_point_3 + 1
    else:
        below_point_3 = below_point_3 + 1         

print('lev AA for different superfamily:total',len(diff_homo), ' above point 0.3', above_point_3 ,'below point 0.3', below_point_3)

exit()


same_homo_df = pd.DataFrame({'AA_seq_idnt':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'AA_seq_idnt':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["AA_seq_idnt"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/seq_idnt_AA_super_safe_cath2.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/seq_idnt_AA_super_safe_cath2.svg', format="svg")
plt.close()
same_homo_not_same_homo_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/AA_same_different_superfamily.csv',index=False)





