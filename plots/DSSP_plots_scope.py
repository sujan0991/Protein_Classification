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

import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, auc,roc_auc_score
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report



df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_DSSP_results.csv')

length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_ForkPoolWorker-1.csv')
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')
length_df = length_df[['id','ss_length']]

print(len(length_df))

family_df = pd.read_csv('/group/bioinf_protstr/article_briefings/data/scope/scope.csv')
family_df = family_df[['index','scop']]
family_df.columns = ['id','scop']

length_df =length_df.merge(family_df, left_on='id', right_on='id', how='inner')

df = pd.merge(df, length_df ,how='left', on="id")
df.columns = ['protein1', 'id', 'DSSP_SS_distance', 'ss_length1', 'scop1']
df = pd.merge(df, length_df ,how='left', on="id")
df=df[~df.scop.isna()]
df.columns = ['protein1','protein2', 'DSSP_SS_distance', 'protein1_length', 'scop1','protein2_length', 'scop2']

df['max_seq_len'] = df[["protein1_length", "protein2_length"]].max(axis=1)
df['levIdent_SS'] = pd.to_numeric(df['DSSP_SS_distance'])/ df['max_seq_len']
df['SS_score_DSSP'] = abs(1-df["levIdent_SS"]).astype(np.float16)


df["scope_superFamily1"] = df['scop1'].str.rpartition('.')[0]
df["scope_superFamily2"] = df['scop2'].str.rpartition('.')[0]

df['is_same_scope_Superfamily'] = np.where(df["scope_superFamily1"] == df["scope_superFamily2"], 1, 0)


df = df[['protein1','protein2','is_same_scope_Superfamily','SS_score_DSSP']]

df['key_id']=['-'.join(sorted(combine)) for combine in zip(df['protein1'], df['protein2'])]




tm_df = pd.read_csv('/group/bioinf_protstr/article_briefings/data/scope/structure/structure.csv')
tm_df["tm_min"] = tm_df[["tm_query", "tm_subject"]].min(axis=1)
tm_df["tm_max"] = tm_df[["tm_query", "tm_subject"]].max(axis=1)
tm_df = tm_df[['query', 'subject', 'tm_min', 'tm_max']]
tm_df['key_id']=['-'.join(sorted(combine)) for combine in zip(tm_df['query'], tm_df['subject'])]


df = df.merge(tm_df,left_on='key_id',right_on='key_id')

df2 = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_results.csv')
df2 = df2[df2.AA_distance != "AA_distance"]

df2['key_id']=['-'.join(sorted(combine)) for combine in zip(df2['id'], df2['query_P'])]


length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_ForkPoolWorker-1.csv')
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')
length_df = length_df[['id','ss_length']]

print(len(length_df))

family_df = pd.read_csv('/group/bioinf_protstr/article_briefings/data/scope/scope.csv')
family_df = family_df[['index','scop']]
family_df.columns = ['id','scop']

length_df =length_df.merge(family_df, left_on='id', right_on='id', how='inner')

###########
##length_df =length_df[['id','ss_length','scop']]

df2 = pd.merge(df2, length_df ,how='left', on="id")

##['id', 'query_P', 'AA_distance', 'SS_distance', 'key_id', 'ss_length_x', 'scop_x', 'ss_length_y', 'scop_y']

df2 = df2[['id', 'query_P','key_id', 'AA_distance', 'SS_distance', 'ss_length', 'scop']]

df2.columns = ['protein1', 'id','key_id', 'AA_distance', 'SS_distance', 'protein1_length', 'scop1']
df2 = pd.merge(df2, length_df ,how='left', on="id")

df2.columns = ['protein1','protein2','key_id','AA_distance', 'SS_distance', 'protein1_length', 'scop1','protein2_length', 'scop2']


df2['max_seq_len'] = df2[["protein1_length", "protein2_length"]].max(axis=1)
df2['levIdent_SS'] = pd.to_numeric(df2['SS_distance'])/ df2['max_seq_len']
df2['SS_score'] = abs(1-df2["levIdent_SS"]).astype(np.float16)
df2['levIdent_AA'] = pd.to_numeric(df2['AA_distance'])/ df2['max_seq_len']
df2['AA_score'] = abs(1-df2["levIdent_AA"]).astype(np.float16)
df2['scope_Family'] = np.where(df2["scop1"] == df2["scop2"], 1, 0)

df2 = df2[['key_id','AA_score','SS_score','scope_Family']]

df = df.merge(df2,left_on='key_id',right_on='key_id')


df = df[['protein1', 'protein2', 'is_same_scope_Superfamily', 'SS_score_DSSP', 'key_id','tm_min', 'tm_max', 'AA_score', 'SS_score', 'scope_Family']]


df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_LEV_TM_results.csv')




#df=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_LEV_TM_results.csv')



###### balance acc

fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['SS_score_DSSP'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2

best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
#balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

plt.axis([0, 1, 0, 1])
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

# plt.axvline(x=thresholds[(thresholds>=0.5)&(thresholds<=0.5)],ymax=round(balanced_accuracy_05[0],2), color="darkgrey",linestyle='--')
# plt.axhline(y=round(balanced_accuracy_05[0],2),xmax=thresholds[(thresholds>=0.5)&(thresholds<=0.5)], color="darkgrey",linestyle='--')
# plt.figtext(0.08, 0.76, round(balanced_accuracy_05[0],2), wrap=True, horizontalalignment='center', fontsize=12) #### ss = 0.74



plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/SS_score_DSSP_superfamily_balanced_accuracy_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/SS_score_DSSP_superfamily_balanced_accuracy_scope.svg', format="svg")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/aa_superfamily_balanced_accuracy_scope.svg', format="svg")

plt.close()





###### latters frequency

# H_alph = ['G', 'H', 'I']
# S_alph = ['E', 'B']
# L_alph = ['S', 'T', 'N']

df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_seq_with_N.csv')

ss_list = df['dssp'].tolist()


h_g_count = 0
h_h_count = 0
h_i_count = 0
s_e_count = 0
s_b_count = 0
l_s_count = 0
l_t_count = 0
l_n_count = 0

for x in ss_list:
    g = x.count('G')
    h = x.count('H')
    i = x.count('I')
    e = x.count('E')
    b = x.count('B')
    s = x.count('S')
    t = x.count('T')
    n = x.count('N')
    h_g_count = h_g_count + g
    h_h_count = h_h_count + h
    h_i_count = h_i_count + i
    s_e_count = s_e_count + e
    s_b_count = s_b_count + b
    l_s_count = l_s_count + s
    l_t_count = l_t_count + t
    l_n_count = l_n_count + n


letters = ['G','H','I','E','B','S','T','N']
h_g_frequency = h_g_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)
h_h_frequency = h_h_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)
h_i_frequency = h_i_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)
s_e_frequency = s_e_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)
s_b_frequency = s_b_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)
l_s_frequency = l_s_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)
l_t_frequency = l_t_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)
l_n_frequency = l_n_count/(h_g_count+h_h_count+h_i_count+s_e_count+s_b_count+l_s_count+l_t_count+l_n_count)

frequency = [h_g_frequency,h_h_frequency,h_i_frequency,s_e_frequency,s_b_frequency,l_s_frequency,l_t_frequency,l_n_frequency]
df = pd.DataFrame({
    'letters': letters,
    'frequency': frequency
})


sns.barplot(data=df, x="letters", y="frequency")
plt.ylabel('Relative Frequency')
plt.xlabel('')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/letter_distribution_SS_DSSP_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/letter_distribution_SS_DSSP_scope.svg', format="svg")
plt.close()

