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
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, auc,roc_auc_score
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report





aa_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/hh_result_processed_all_safe_cath_after_ck1-2-3.csv',usecols =['key_id','min_I','max_I'])

blast_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/BLAST_CATH/blast.csv')
blast_df['blast_identity'] = (blast_df['%id']/100)* blast_df['qcov']
blast_df = blast_df[['query','subject','blast_identity','qcov','E_value']]
blast_df.columns = ['id','subject','blast_identity','qcov','E_value']



###################

lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_TM_results_all_safe_cath_after_ck1-2-3.csv')

lev_tm_df = lev_tm_df.merge(aa_df,left_on='key_id',right_on='key_id')

blast_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/safe_cath_after_ck1-2-3_blast_result_high_e_value.csv')
blast_df.columns = ['query','subject','%id','alignment_length','mismatches','gap_openings','query_start','query_end','subject_start','subject_end','E_value','bit_score','qcov']
blast_df = blast_df[['query','subject','%id','qcov','E_value']]
blast_df['key_id']=['-'.join(sorted(combine)) for combine in zip(blast_df['query'], blast_df['subject'])]
blast_df.drop_duplicates(subset='key_id', keep="last", inplace=True)

lev_tm_df = pd.merge(lev_tm_df, blast_df ,how='left', on="key_id")


##lev_tm_df.isna().any()
lev_tm_df=lev_tm_df[~lev_tm_df.blast_identity_min.isna()]



lev_tm_df["tm_min_class"] = np.where(lev_tm_df["TM_min"] >=0.5, 1, 0)
lev_tm_df["seq_idnt_class_SS"] = np.where(lev_tm_df['SS_score'] >=0.5, 1, 0)
lev_tm_df["seq_idnt_class_AA"] = np.where(lev_tm_df['AA_score'] >=0.3, 1, 0)
#lev_tm_df["cath_superFamily"] = np.where(lev_tm_df["Homol1"] == lev_tm_df["Homol2"], 1, 0)


#### although i ploted this three polt in merge classes, i am doing it again

##  SS

same_homo = lev_tm_df.query('cath_superFamily == 1') ['SS_score']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['SS_score']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.7:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 0.5', above_point_5 ,'below point 0.5', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.7:
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
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_super_have_blast_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_super_have_blast_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


## AA

same_homo = lev_tm_df.query('cath_superFamily == 1') ['AA_score']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['AA_score']

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
 
same_homo_df = pd.DataFrame({'AA_score':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'AA_score':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["AA_score"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/AA_score_super_have_blast_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/AA_score_super_have_blast_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


### TM

same_homo = lev_tm_df.query('cath_superFamily == 1') ['TM_max']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['TM_max']

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
 
same_homo_df = pd.DataFrame({'TM_max':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'TM_max':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["TM_max"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_all_data_score_super_have_blast_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_all_data_score_super_have_blast_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()

###########################################################################

### min idendity
same_homo = lev_tm_df.query('cath_superFamily == 1') ['min_I']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['min_I']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 50:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 50', above_point_5 ,'below point 50', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 50:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 50', above_point_5 ,'below point 50', below_point_5)
 


same_homo_df = pd.DataFrame({'min_I':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'min_I':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["min_I"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/Identity_min_score_super_have_blast_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/Identity_min_score_super_have_blast_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()

### max idendity
same_homo = lev_tm_df.query('cath_superFamily == 1') ['max_I']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['max_I']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 30:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 30', above_point_5 ,'below point 30', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 30:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 30', above_point_5 ,'below point 30', below_point_5)
 


same_homo_df = pd.DataFrame({'max_I':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'max_I':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["max_I"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/Identity_max_score_super_have_blast_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/Identity_max_score_super_have_blast_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()

## blast



### min

same_homo = lev_tm_df.query('cath_superFamily == 1') ['blast_identity_min']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['blast_identity_min']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 30:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 30', above_point_5 ,'below point 30', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 30:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 30', above_point_5 ,'below point 30', below_point_5)
 

same_homo_df = pd.DataFrame({'blast_identity_min':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity_min':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity_min"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_min_all_data_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_min_all_data_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


### max

same_homo = lev_tm_df.query('cath_superFamily == 1') ['blast_identity_max']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['blast_identity_max']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 30:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 30', above_point_5 ,'below point 30', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 30:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 30', above_point_5 ,'below point 30', below_point_5)
 

same_homo_df = pd.DataFrame({'blast_identity_max':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity_max':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity_max"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_max_all_data_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_max_all_data_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


####### 

###### blast_identity_min*qcov_min

lev_tm_df['blast_identity_min*qcov_min'] = lev_tm_df['blast_identity_min']* lev_tm_df['qcov_min']


same_homo = lev_tm_df.query('cath_superFamily == 1') ['blast_identity_min*qcov_min']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['blast_identity_min*qcov_min']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 2000:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 2000', above_point_5 ,'below point 2000', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 2000:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 2000', above_point_5 ,'below point 2000', below_point_5)
 


same_homo_df = pd.DataFrame({'blast_identity_min*qcov_min':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity_min*qcov_min':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df['blast_identity_min*qcov_min'], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_min_cov_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_min_cov_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


###### blast_identity_max*qcov_max

lev_tm_df['blast_identity_max*qcov_max'] = lev_tm_df['blast_identity_max']* lev_tm_df['qcov_max']

same_homo = lev_tm_df.query('cath_superFamily == 1') ['blast_identity_max*qcov_max']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['blast_identity_max*qcov_max']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 5000:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 5000', above_point_5 ,'below point 5000', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 5000:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 5000', above_point_5 ,'below point 5000', below_point_5)
 



same_homo_df = pd.DataFrame({'blast_identity_max*qcov_max':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity_max*qcov_max':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity_max*qcov_max"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_max_cov_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_max_cov_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()



same_homo = lev_tm_df.query('cath_superFamily == 1') ['blast_identity_max*qcov']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['blast_identity_max*qcov']
same_homo_df = pd.DataFrame({'blast_identity_max*qcov':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity_max*qcov':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity_max*qcov"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/blast_identity_max_cov_2_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/blast_identity_max_cov_2_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


###




lev_tm_df['E_value_min'].isnull().values.any()
lev_tm_df['E_value_min'].isnull().sum()

same_homo = lev_tm_df.query('cath_superFamily == 1') ['E_value_min']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['E_value_min']

same_homo_df = pd.DataFrame({'E_value_min':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'E_value_min':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["E_value_min"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_E_value_min_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_E_value_min_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()




#################

##lev_tm_df["tm_min_round"] = lev_tm_df["TM_min"].round(decimals = 1)
seq_dist_SS = lev_tm_df['SS_score'].to_list()
tm_score_min = lev_tm_df["TM_min"].to_list()

tm_0_1 = 0
tm_1_2 =0
tm_2_3 = 0
tm_3_4 = 0
tm_4_5 =0
tm_5_6 = 0
tm_6_7 = 0
tm_7_8 =0
tm_8_9 = 0
tm_9_1 = 0



for i,value in enumerate(tm_score_min):
    if tm_score_min[i] >= 0.0 and tm_score_min[i] < 0.1:
        tm_0_1 = tm_0_1 + 1
        tm_score_min[i] = 0.05
    elif tm_score_min[i] >= 0.1 and tm_score_min[i] < 0.2:
        tm_1_2 = tm_1_2 + 1
        tm_score_min[i] = 0.15
    elif tm_score_min[i] >= 0.2 and tm_score_min[i] < 0.3:
        tm_2_3 = tm_2_3 + 1
        tm_score_min[i] = 0.25
    elif tm_score_min[i] >= 0.3 and tm_score_min[i] < 0.4:
        tm_3_4 = tm_3_4 + 1
        tm_score_min[i] = 0.35
    elif tm_score_min[i] >= 0.4 and tm_score_min[i] < 0.5:
        tm_4_5 = tm_4_5 + 1
        tm_score_min[i] = 0.45
    elif tm_score_min[i] >= 0.5 and tm_score_min[i] < 0.6:
        tm_5_6 = tm_5_6 + 1
        tm_score_min[i] = 0.55
    elif tm_score_min[i] >= 0.6 and tm_score_min[i] < 0.7:
        tm_6_7 = tm_6_7 + 1
        tm_score_min[i] = 0.65
    elif tm_score_min[i] >= 0.7 and tm_score_min[i] < 0.8:
        tm_7_8 = tm_7_8 + 1
        tm_score_min[i] = 0.75
    elif tm_score_min[i] >= 0.8 and tm_score_min[i] < 0.9:
        tm_8_9 = tm_8_9 + 1
        tm_score_min[i] = 0.85  
    else :                               
        ##tm_score_min[i] >= 0.9 and tm_score_min[i] < 1.0:
        tm_9_1 = tm_9_1 + 1
        tm_score_min[i] = 0.95

print("tm-min 0-1",tm_0_1,tm_1_2,tm_2_3,tm_3_4,tm_4_5,tm_5_6,tm_6_7,tm_7_8,tm_8_9,tm_9_1) 

temp_dict = {'tm_score_min':tm_score_min, 'SS_score':seq_dist_SS}
temp_df = pd.DataFrame(temp_dict)
#print("temp_df........len",len(temp_df))

plt.figure(figsize=(20,8))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
sns.violinplot(data=temp_df, x='tm_score_min', y='SS_score')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_SS_scale_0_1_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_SS_scale_0_1_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()



############


lev_tm_df['max_identities_class_03']=(lev_tm_df.max_I>=30) &(~lev_tm_df.max_I.isna())
lev_tm_df['max_identities_class_02']=(lev_tm_df.max_I>=20) &(~lev_tm_df.max_I.isna())
tn, fp, fn, tp = confusion_matrix(lev_tm_df.cath_superFamily , lev_tm_df['max_identities_class_03']).ravel()
print('max_identities_class_03','TNR, FPR,  FNR, TPR', (tn/(tn+fp)), (fp/(fp+tn)), (fn/(tp+fn)), (tp/(tp+fn)))

lev_tm_df['min_identities_class_03']=(lev_tm_df.min_I>=30) &(~lev_tm_df.min_I.isna())
lev_tm_df['min_identities_class_02']=(lev_tm_df.min_I>=20) &(~lev_tm_df.min_I.isna())
tn, fp, fn, tp = confusion_matrix(lev_tm_df.cath_superFamily , lev_tm_df['min_identities_class_03']).ravel()
print('min_identities_class_03','TNR, FPR,  FNR, TPR', (tn/(tn+fp)), (fp/(fp+tn)), (fn/(tp+fn)), (tp/(tp+fn)))




##


lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_lev_tm_HH_blast_blast_ss_results.csv',usecols=['cath_superFamily','TM_min','SS_score'])


########## AUC

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['TM_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...tm.max",roc_auc)

plt.plot(fpr,tpr,color="blue",label="TM_max")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['SS_score'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))

#create ROC curve
plt.plot(fpr,tpr,color="green",label="SS")


# fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['min_I'], pos_label=1)
# print("fpr, tpr len",len(fpr),len(tpr))
# roc_auc = auc(fpr, tpr)
# print("roc_auc...min_I.",roc_auc)

# plt.plot(fpr,tpr,color="orange",label="min_I")


plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_max_lev_min_I_roc_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/auc/svg/tm_max_lev_roc_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


##########

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['E_value_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...E_value_min.",roc_auc)

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['E_value_max'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...E_value_max.",roc_auc)


####### balanced_accuracy 

## E_value_min

#fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['SS_score'], pos_label=1)
#fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['AA_score'], pos_label=1)
#fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['TM_min'], pos_label=1)
fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['max_I'], pos_label=1)
#fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['E_value_max'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2

##### for table 

# df_from_arr = pd.DataFrame(data=[fpr, tpr,thresholds,tnr,fnr,balanced_accuracy]).T
# df_from_arr.columns = ['fpr','tpr','thresholds','tnr','fnr','balanced_accuracy']
# df_from_arr.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/balanced_accuracy_table_childrenFilter_False_v3_onlyHomo.csv',index=False)

############

best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

#plt.axis([0, 1, 0, 1])
plt.plot(thresholds,tpr, color="blue",label="tpr" )
#plt.plot(thresholds,fpr, color="purple",label="fpr" )
plt.plot(thresholds,tnr, color="green", label="tnr")
#plt.plot(thresholds,fnr, color="orange",label="fnr" )
plt.plot(thresholds,balanced_accuracy, color="black",label="ba")
plt.legend(loc="lower right")

#plt.axvline(x=best_thre[0],ymax=max(balanced_accuracy), color="darkgrey",linestyle='--')
plt.figtext(0.41, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12) ## for aa x=0.30,ss=0.55
#plt.axhline(y=round(max(balanced_accuracy),2),xmax=round(best_thre[0],2), color="darkgrey",linestyle='--')
plt.figtext(0.08, 0.73, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)## for aa y=0.59,ss= 0.78


# plt.savefig('tm_min_superfamily_balanced_accuracy_test_set_2.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/max_I_superfamily_balanced_accuracy_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/max_I_superfamily_balanced_accuracy_safe_cath_after_ck1-2-3.svg', format="svg")

plt.close()



### tm == 0.5

tm_05 = lev_tm_df.query('TM_min == 0.5')



## data distribution of tm,ss,aa

fig = sns.kdeplot(data=lev_tm_df, x="TM_min", shade=True, color="g")
fig = sns.kdeplot(data=lev_tm_df, x="SS_seq_idnt", shade=True, color="b")
fig = sns.kdeplot(data=lev_tm_df, x="AA_seq_idnt", shade=True, color="m")

plt.legend(loc="upper right")
plt.legend(labels=["TM","SS","AA"])
plt.xlim(0,1)
plt.xlabel('Score')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/TM_SS_AA_distribution_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/TM_SS_AA_distribution_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()




## data distribution of  same and different superfamily

# same_sf = lev_tm_df.query('Homol1 == Homol2')

# fig = sns.kdeplot(data=same_sf, x="TM_min", shade=True, color="g")
# fig = sns.kdeplot(data=same_sf, x="SS_seq_idnt", shade=True, color="b")
# fig = sns.kdeplot(data=same_sf, x="AA_seq_idnt", shade=True, color="m")

# plt.legend(loc="upper right")
# plt.legend(labels=["TM","SS","AA"])
# plt.xlim(0,1)
# plt.xlabel('Score')
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/same_uperfamily_distribution_TM_SS_AA_safe_cath_after_ck1-2-3.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/same_uperfamily_distribution_TM_SS_AA_csafe_cath_after_ck1-2-3.svg', format="svg")
# plt.close()



# diff_sf = lev_tm_df.query('Homol1 != Homol2')

# fig = sns.kdeplot(data=diff_sf, x="TM_min", shade=True, color="g")
# fig = sns.kdeplot(data=diff_sf, x="SS_seq_idnt", shade=True, color="b")
# fig = sns.kdeplot(data=diff_sf, x="AA_seq_idnt", shade=True, color="m")

# plt.legend(loc="upper right")
# plt.legend(labels=["TM","SS","AA"])
# plt.xlim(0,1)
# plt.xlabel('Score')
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/different_uperfamily_distribution_TM_SS_AA_safe_cath_after_ck1-2-3.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/different_uperfamily_distribution_TM_SS_AA_safe_cath_after_ck1-2-3.svg', format="svg")
# plt.close()





####### TM vs SS

# seq_dist_SS = lev_tm_df['seq_idnt'].to_list()
# seq_dist_AA = lev_tm_df['AA_seq_idnt'].to_list()
# tm_score_min = lev_tm_df["TM_min"].to_list()
# seq_dist_SS = np.nan_to_num(seq_dist_SS)
# tm_score_min = np.nan_to_num(tm_score_min)
# seq_dist_AA = np.nan_to_num(seq_dist_AA)


# temp_dict = {'tm_score':tm_score_min, 'seq_idnt_SS':seq_dist_SS}
# temp_df = pd.DataFrame(temp_dict)


# g=sns.kdeplot(data=temp_df, x="tm_score", y="seq_idnt_SS",fill=True)
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/tm_seq_idnt_SS_den_safe_safe_cath_after_ck1-2-3.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/tm_seq_idnt_SS_den_safe_safe_cath_after_ck1-2-3.svg', format="svg")
# plt.close()



# temp_dict = {'tm_score':tm_score_min, 'seq_idnt_AA':seq_dist_AA}
# temp_df = pd.DataFrame(temp_dict)


# g=sns.kdeplot(data=temp_df, x="tm_score", y="seq_idnt_AA",fill=True)
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/tm_seq_idnt_AA_den_safe_cath2.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/tm_seq_idnt_AA_den_safe_cath2.svg', format="svg")
# plt.close()






### all info

# df_merged=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/merged_Cath_list_cleaned.csv')

# ## 31893
# test_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_ForkPoolWorker-2.csv')
# test_df.drop_duplicates(subset='id', keep="last", inplace=True)
# test_df=test_df[~test_df.SS_seq.isna()]
# test_df=test_df[~test_df.AA_seq.isna()]
# test_df=test_df[~test_df.id.isna()]
# test_df.reset_index(inplace=True)


length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv', usecols=['id','AA_seq','SS_seq'])
length_df.drop_duplicates(subset='id', keep="last", inplace=True)
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)

df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')


df_24k = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/domain_id_list_24k.csv')
list_24k = df_24k['domain_id'].tolist()

length_df = length_df[length_df['id'].isin(list_24k)]
length_df.isna().any()

## ss letter distribution

ss_list = length_df['SS_seq'].tolist()



h_count = 0
s_count = 0
l_count = 0

for x in ss_list:
    h = x.count('H')
    s = x.count('S')
    l = x.count('L')
    h_count = h_count + h
    s_count = s_count + s
    l_count = l_count + l

print("h_count,s_count l_count :",h_count,s_count,l_count)    

letters = ['H','S','L']
h_frequency = h_count/(h_count+s_count+l_count)
s_frequency = s_count/(h_count+s_count+l_count)
l_frequency = l_count/(h_count+s_count+l_count)
frequency = [h_frequency,s_frequency,l_frequency]
df = pd.DataFrame({
    'letters': letters,
    'frequency': frequency
})


sns.barplot(data=df, x="letters", y="frequency")
plt.ylabel('Relative Frequency')
plt.xlabel('')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/letter_distribution_SS_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/letter_distribution_SS_after_ck1-2-3.svg', format="svg")
plt.close()
      
    


# test_df =test_df[['id']] ## 27842
# test_df = test_df['id'].tolist()



# df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-domain-description-v4_3_0.csv')
# print(len(df_family.Homol.unique()))

# print(len(pd.concat([lev_tm_df['Homol1'], lev_tm_df['Homol2']]).unique()))


length_df = length_df[length_df.ss_length > 50] ## it seems like some pdb have problem with the selection range, 
                                                  ##thats why the length check in SuperFamily_20k didn't work for those## example : 1ahsA00

sns.displot(length_df, x="ss_length", kde=True)
plt.xlabel('Number of Residues')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/length_distribution_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/length_distribution_after_ck1-2-3.svg', format="svg")
plt.close()





# df_merged = df_merged[df_merged.id.isin(test_df)]
# sns.displot(df_merged, x="length", kde=True)
# plt.xlabel('Number of Residues')
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/length_distribution_batch_0_safe_cath2.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/length_distribution_batch_0_safe_cath2.svg', format="svg")
# plt.close()


# ## length difference plot
# # aa=([np.log10(i) for i in helix.values()])  

# lev_tm_df['length_diff'] = abs(lev_tm_df['domain1_length'] - lev_tm_df['domain2_length'])

# sns.displot(lev_tm_df, x="length_diff", kde=True)
# ##plt.figure(figsize=(12,10))
# plt.yscale("log") 
# plt.xlabel('Length Difference')
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/length_diff_distribution_batch_0_safe_cath2.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/length_diff_distribution_batch_0_safe_cath2.svg', format="svg")
# plt.close()




#########

length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv', usecols=['id','AA_seq','SS_seq'])
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')
#length_df = length_df[['id','ss_length']]

print(len(length_df))

df_24k = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/domain_id_list_24k.csv')
list_24k = df_24k['domain_id'].tolist()

length_df = length_df[length_df['id'].isin(list_24k)]

df2 = length_df.query('AA_seq.str.contains("\?")')
df2_ids = df2['id'].tolist()

length_df = length_df[~length_df['id'].isin(df2_ids)]



df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')


lev_tm_df.columns = ['id', 'domain2', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max']
lev_tm_df = pd.merge(lev_tm_df, length_df ,how='left', on="id")
lev_tm_df.columns = ['domain1', 'id', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max', 'Class1', 'Arch1', 'Topol1', 'Homol1']
lev_tm_df = pd.merge(lev_tm_df, length_df ,how='left', on="id")
lev_tm_df.columns = ['domain1', 'domain1', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'Class2', 'Arch2', 'Topol2', 'Homol2']

##lev_tm_df.isna().any()

same_homo = lev_tm_df.query('Homol1 == Homol2')
diff_homo = lev_tm_df.query('Homol1 != Homol2')

print('same',len(same_homo),'different',len(diff_homo))





###### class distribution

length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv', usecols=['id','AA_seq','SS_seq'])
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')
length_df = length_df[['id','ss_length']]

print(len(length_df))

df_24k = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/domain_id_list_24k.csv')
list_24k = df_24k['domain_id'].tolist()

length_df = length_df[length_df['id'].isin(list_24k)]



df2 = length_df.query('AA_seq.str.contains("\?")')
df2_ids = df2['id'].tolist()

length_df = length_df[~length_df['id'].isin(df2_ids)]

df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')


length_df = length_df[['id','Class']]

cl_count = length_df.Class.value_counts()

df = pd.DataFrame({'class':cl_count.index, 'count':cl_count.values})


plt.figure(figsize=(10,8))
sns.barplot(data=df, x="class", y="count")
plt.ylabel('Count')
plt.xlabel('')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/class_distribution_cath_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/class_distribution_cath_ck1-2-3.svg', format="svg")
plt.close()




## resi len distribution


df_cath=pd.read_csv("/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-selection-v4_3_0.csv")
df_cath['length']= df_cath['selection'].str.split(' resi ',expand=True)[1]
df_cath = df_cath[['id','length']]
##df_cath.loc[df_cath['comma_count'].idxmax()]

df = pd.concat([df_cath['id'], df_cath['length'].str.split(',', expand = True)], axis = 1)
df_cath = pd.merge(df_cath, df ,how='left', on="id")
df_cath = df_cath.replace(['NaN', 'None', ''], float('nan'))
df_cath = df_cath.fillna('0-0')

df_cath['sel0'] = abs((df_cath[0].str.split('-',expand=True)[0]).astype(int) - (df_cath[0].str.split('-',expand=True)[1]).astype(int)) + 1
df_cath['sel1'] = abs((df_cath[1].str.split('-',expand=True)[0]).astype(int) - (df_cath[1].str.split('-',expand=True)[1]).astype(int)) + 1
df_cath['sel2'] = abs((df_cath[2].str.split('-',expand=True)[0]).astype(int) - (df_cath[2].str.split('-',expand=True)[1]).astype(int)) + 1
df_cath['sel3'] = abs((df_cath[3].str.split('-',expand=True)[0]).astype(int) - (df_cath[3].str.split('-',expand=True)[1]).astype(int)) + 1
df_cath['sel4'] = abs((df_cath[4].str.split('-',expand=True)[0]).astype(int) - (df_cath[4].str.split('-',expand=True)[1]).astype(int)) + 1
df_cath['sel5'] = abs((df_cath[5].str.split('-',expand=True)[0]).astype(int) - (df_cath[5].str.split('-',expand=True)[1]).astype(int)) + 1

df_cath = df_cath.replace(1, 0)

column_names = ['sel0', 'sel1', 'sel2', 'sel3', 'sel4', 'sel5']
df_cath['total_length']= df_cath[column_names].sum(axis=1)

df_cath.loc[df_cath['total_length'].idxmin()]
df_cath.loc[df_cath['total_length'].idxmax()]
df_cath.loc[:, 'total_length'].mean()


sns.displot(df_cath, x="total_length", kde=True)
plt.xlabel('Number of Residues')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_all_cath.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_all_cath.svg', format="svg")
plt.close()


############ test

# df_cath= pd.DataFrame({'id': ['9xiaA00','9xiaA01','103lA00','9ximC00','9ximA00'],
#                   'length': ['1-163','1-163,164-170','1-160,164-170,171-175','1-163,164-170,171-175,176-180','1-163,164-170,171-175,176-180,181-190,191-200']
#                 })

# df = pd.concat([df_cath['id'], df_cath['length'].str.split(',', expand = True)], axis = 1)
# df_cath = pd.merge(df_cath, df ,how='left', on="id")
# df_cath = df_cath.replace(['NaN', 'None', ''], float('nan'))
# df_cath = df_cath.fillna('0-0')

# df_cath['sel0'] = abs((df_cath[0].str.split('-',expand=True)[0]).astype(int) - (df_cath[0].str.split('-',expand=True)[1]).astype(int)) + 1
# df_cath['sel1'] = abs((df_cath[1].str.split('-',expand=True)[0]).astype(int) - (df_cath[1].str.split('-',expand=True)[1]).astype(int)) + 1
# df_cath['sel2'] = abs((df_cath[2].str.split('-',expand=True)[0]).astype(int) - (df_cath[2].str.split('-',expand=True)[1]).astype(int)) + 1
# df_cath['sel3'] = abs((df_cath[3].str.split('-',expand=True)[0]).astype(int) - (df_cath[3].str.split('-',expand=True)[1]).astype(int)) + 1
# df_cath['sel4'] = abs((df_cath[4].str.split('-',expand=True)[0]).astype(int) - (df_cath[4].str.split('-',expand=True)[1]).astype(int)) + 1
# df_cath['sel5'] = abs((df_cath[5].str.split('-',expand=True)[0]).astype(int) - (df_cath[5].str.split('-',expand=True)[1]).astype(int)) + 1

# df_cath = df_cath.replace(1, 0)

# column_names = ['sel0', 'sel1', 'sel2', 'sel3', 'sel4', 'sel5']
# df_cath['total_length']= df_cath[column_names].sum(axis=1)




####  exclude anything that does not contain any of these patterns: HLH, HLS, SLS, SLH
## And provide the above table.

## bash command
### perl -nE 'chomp; $x=0; $x=1 if m/.*(H.*L.*H)|(H.*S.*L)|(S.*L.*S)|(S.*L.*H)|(L.*H.*S)|(L.*S.*H).*/; say $_.",".$x' /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/cath_24k_SS_sequences.csv  | grep ',0' 


df = pd.read_csv('/home/bioinf/mdho200b/bad_cath.txt')

bad_id_list = df['id'].tolist()

lev_tm_df = lev_tm_df[~lev_tm_df['domain1'].isin(bad_id_list)]
lev_tm_df = lev_tm_df[~lev_tm_df['domain2'].isin(bad_id_list)]

lev_tm_df = lev_tm_df[['domain1', 'domain2', 'cath_superFamily', 'SS_score', 'AA_score', 'key_id', 'TM_min', 'TM_max', 'min_I', 'max_I', 'blast_identity_min', 
                       'blast_identity_max', 'qcov_min', 'qcov_max', 'E_value_min', 'E_value_max']]


###lev_tm_df = lev_tm_df.merge(df,left_on='key_id',right_on='key_id')  ## this df is dssp df



########## AUC

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['TM_max'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...tm.max",roc_auc)


fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['SS_score'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))



fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['AA_score'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc....aa",roc_auc)


fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['min_I'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...min_I.",roc_auc)



fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['max_I'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...max_I.",roc_auc)


lev_tm_df['blast_identity_min'] = lev_tm_df['blast_identity_min'].fillna(0)

lev_tm_df['blast_identity_min'].isnull().values.any()

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['blast_identity_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity_min.",roc_auc)



lev_tm_df['blast_identity_max'] = lev_tm_df['blast_identity_max'].fillna(0)
lev_tm_df['blast_identity_max'].isnull().values.any()

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['blast_identity_max'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity_max.",roc_auc)

lev_tm_df['E_value_min'] = lev_tm_df['E_value_min'].fillna(1)
lev_tm_df['E_value_max'] = lev_tm_df['E_value_min'].fillna(1)
lev_tm_df['E_value_min'] = lev_tm_df['E_value_min'] + 0.000000001
lev_tm_df['E_value_max'] = lev_tm_df['E_value_max'] + 0.000000001

lev_tm_df['E_value_min'].isnull().values.any()
lev_tm_df['E_value_max'].isnull().values.any()

lev_tm_df['E_value_min'] = np.log10(lev_tm_df['E_value_min'])
lev_tm_df['E_value_max'] = np.log10(lev_tm_df['E_value_max'])

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['E_value_max'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...E_value_max.",roc_auc)


fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['SS_score_DSSP'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...SS_score_DSSP.",roc_auc)



### balance accuracy

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['TM_max'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2


best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

plt.axis([0, 1, 0, 1])
plt.plot(thresholds,tpr, color="blue",label="tpr" )
plt.plot(thresholds,fpr, color="purple",label="fpr" )
plt.plot(thresholds,tnr, color="green", label="tnr")
plt.plot(thresholds,fnr, color="orange",label="fnr" )
plt.plot(thresholds,balanced_accuracy, color="black",label="ba")
plt.legend(loc="lower right")

plt.figtext(0.41, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12) ## for aa x=0.30,ss=0.55
plt.figtext(0.08, 0.73, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)## for aa y=0.59,ss= 0.78

plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_superfamily_balanced_accuracy_safe_cath_after_removing_HLS.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_superfamily_balanced_accuracy_safe_cath_after_removing_HLS.svg', format="svg")

plt.close()






#######



df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv')
#df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_random_1000_speed_test_ForkPoolWorker-12.csv')

df.drop_duplicates(subset='id', keep="last", inplace=True)
df=df[~df.SS_seq.isna()]
df=df[~df.AA_seq.isna()]
df=df[~df.id.isna()]
df.reset_index(inplace=True)

print("....",len(df))

df_24k = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/domain_id_list_24k.csv')
list_24k = df_24k['domain_id'].tolist()

df = df[df['id'].isin(list_24k)]



df2 = df.query('AA_seq.str.contains("\?")')
df2_ids = df2['id'].tolist()

df = df[~df['id'].isin(df2_ids)]



df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
df =df.merge(df_family, left_on='id', right_on='id', how='inner')
df = df[['id', 'Homol']]

fm = df.Homol.value_counts()
count_df = pd.DataFrame({'Homol':fm.index, 'count':fm.values})

df = df.merge(count_df,left_on='Homol',right_on='Homol')


df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/superfamily_count_cath_after_ck1-2-3.csv')