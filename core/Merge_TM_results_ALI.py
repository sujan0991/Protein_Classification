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


####### merge

# path= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_TM_childrenFilter_False_v3/' # use your path
# all_files = glob.glob((path+'*.csv'))
# li = []
# for filename in all_files:
#     df = pd.read_csv(filename, index_col=None, header=0)
#     li.append(df)

# df_merged = pd.concat(li, axis=0, ignore_index=True)
# df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_results_all_childrenFilter_False_v3.csv')

# print('len(df_merged)',len(df_merged))

########

df_merged = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_results_all_childrenFilter_False_v3.csv')
df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-domain-description-v4_3_0.csv')
df_family =df_family[['id','Homol']]

df_merged = pd.merge(df_merged, df_family ,how='left', on="id")
df_merged = df_merged[['id', 'query_P', 'TM_min', 'TM_max', 'Homol']]
df_merged.columns = ['domain1','id','TM_min', 'TM_max','Homol']
df_merged = pd.merge(df_merged, df_family ,how='left', on="id")
#print(df_merged.isnull().any())
df_merged =  df_merged[['domain1','id','TM_min', 'TM_max', 'Homol_x', 'Homol_y']]
df_merged.columns = ['domain1','domain2','TM_min', 'TM_max', 'Homol1', 'Homol2']

df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_results_all_with_superfamily_childrenFilter_False_v3.csv',index=False)


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from scipy.stats import pearsonr

same_homo = df_merged.query('Homol1 == Homol2') ['TM_min']
diff_homo = df_merged.query('Homol1 != Homol2') ['TM_min']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('TM_min for same superfamily:total',len(same_homo), 'above point 0.5', above_point_5 ,'below point 0.5', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('TM_min for different superfamily:total',len(diff_homo), ' above point 0.5', above_point_5 ,'below point 0.5', below_point_5)
 
same_homo_df = pd.DataFrame({'TM_min':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'TM_min':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["TM_min"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_min_super_safe_cath2.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_min_super_safe_cath2.svg', format="svg")
plt.close()
same_homo_not_same_homo_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_min_same_different_superfamily.csv',index=False)

#####


## TMmax

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from scipy.stats import pearsonr

same_homo = df_merged.query('Homol1 == Homol2') ['TM_max']
diff_homo = df_merged.query('Homol1 != Homol2') ['TM_max']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('TM_max for same superfamily:total',len(same_homo), 'above point 0.5', above_point_5 ,'below point 0.5', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('TM_max for different superfamily:total',len(diff_homo), ' above point 0.5', above_point_5 ,'below point 0.5', below_point_5)
 
same_homo_df = pd.DataFrame({'TM_max':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'TM_max':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["TM_max"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_max_super_safe_cath2.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_max_super_safe_cath2.svg', format="svg")
plt.close()
same_homo_not_same_homo_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_max_same_different_superfamily.csv',index=False)

df_merged.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/TM_results_all_with_superfamily_childrenFilter_False_v3_onlyHomo.csv',index=False)
