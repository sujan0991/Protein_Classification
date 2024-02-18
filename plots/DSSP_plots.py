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


df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_DSSP_LEV_results.csv')

length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv', usecols=['id','AA_seq','SS_seq'])
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')
length_df = length_df[['id','ss_length']]

print(len(length_df))

df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')


##length_df =length_df[['id','ss_length','Homol']]

df = pd.merge(df, length_df ,how='left', on="id")
df = df[['id', 'query_P', 'DSSP_SS_distance','ss_length', 'Homol']]
df.columns = ['domain1','id','DSSP_SS_distance', 'domain1_length','Homol']
df = pd.merge(df, length_df ,how='left', on="id")
df = df[['domain1', 'id', 'DSSP_SS_distance', 'domain1_length', 'Homol_x', 'ss_length', 'Homol_y']]
df.columns = ['domain1','domain2','DSSP_SS_distance','domain1_length', 'Homol1','domain2_length', 'Homol2']

df=df[~df.domain1_length.isna()]
df=df[~df.Homol1.isna()]
df=df[~df.domain2_length.isna()]
df=df[~df.Homol2.isna()]




# same_homo = df.query('Homol1 == Homol2') 
# diff_homo = df.query('Homol1 != Homol2')
# print(len(same_homo),len(diff_homo))


df['key_id']=['-'.join(sorted(combine)) for combine in zip(df['domain1'], df['domain2'])]


df['max_seq_len'] = df[["domain1_length", "domain2_length"]].max(axis=1)
df['levIdent_SS'] = pd.to_numeric(df['DSSP_SS_distance'])/ df['max_seq_len']
df['SS_score'] = abs(1-df["levIdent_SS"]).astype(np.float16)
df['cath_superFamily'] = np.where(df["Homol1"] == df["Homol2"], 1, 0)
    #df = df[['SS_score','cath_superFamily']]
df = df[['domain1','domain2', 'cath_superFamily','SS_score','key_id']]
df.columns = ['domain1','domain2', 'cath_superFamily','SS_score_DSSP','key_id']


lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_TM_results_all_safe_cath_after_ck1-2-3.csv')
lev_tm_df = lev_tm_df[['SS_score', 'AA_score', 'key_id', 'TM_min', 'TM_max']]

df = df.merge(lev_tm_df,on='key_id',how='left')




same_homo = df.query('cath_superFamily == 1') ['SS_score_DSSP']
diff_homo = df.query('cath_superFamily == 0') ['SS_score_DSSP']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.7:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 0.7', above_point_5 ,'below point 0.7', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.7:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 0.7', above_point_5 ,'below point 0.7', below_point_5)


same_homo_df = pd.DataFrame({'SS_score':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'SS_score':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot(y=same_homo_not_same_homo_df["SS_score"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_DSSP_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_DSSP_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


### AA 

same_homo = df.query('cath_superFamily == 1') ['AA_score']
diff_homo = df.query('cath_superFamily == 0') ['AA_score']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.3:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 0.3', above_point_5 ,'below point 0.3', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.3:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 0.3', above_point_5 ,'below point 0.3', below_point_5)


same_homo_df = pd.DataFrame({'AA_score':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'AA_score':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot(y=same_homo_not_same_homo_df["AA_score"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/AA_score_DSSP_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/AA_score_DSSP_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


### tm 

same_homo = df.query('cath_superFamily == 1') ['TM_max']
diff_homo = df.query('cath_superFamily == 0') ['TM_max']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.3:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 0.3', above_point_5 ,'below point 0.3', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.3:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 0.3', above_point_5 ,'below point 0.3', below_point_5)


same_homo_df = pd.DataFrame({'TM_max':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'TM_max':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot(y=same_homo_not_same_homo_df["TM_max"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_score_DSSP_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_score_DSSP_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()




######## AUC

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['SS_score_DSSP'], pos_label=1)
print("fpr, tpr len",fpr,tpr)
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))


plt.plot(fpr,tpr,color="green",label="SS-DSSP")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['AA_score'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc....aa",roc_auc)

plt.plot(fpr,tpr,color="orange",label="AA")
#plt.plot([], [], ' ', label= round(roc_auc, 2))

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['TM_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...tm.min",roc_auc)

plt.plot(fpr,tpr,color="blue",label="TM_min")


fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['TM_max'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...tm.max",roc_auc)

plt.plot(fpr,tpr,color="cyan",label="TM_max")


plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_max_min_lev_roc_DSSP_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_max_min_lev_roc_DSSP_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()



# plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

# plt.legend(loc="lower right")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/ss_roc_DSSP_safe_cath_after_ck1-2-3.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/ss_roc_DSSP_safe_cath_after_ck1-2-3.svg', format="svg")
# plt.close()



############    balance accuracy

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['SS_score_DSSP'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2

best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

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



plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/ss_DSSP_superfamily_balanced_accuracy_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/ss_DSSP_superfamily_balanced_accuracy_safe_cath_after_ck1-2-3.svg', format="svg")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/seq_idnt_AA_superfamily_balanced_accuracy_safe_cath2.svg', format="svg")

plt.close()




####


seq_dist_SS = df['SS_score_DSSP'].to_list()
tm_score_min = df["TM_min"].to_list()

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

temp_dict = {'tm_score_min':tm_score_min, 'SS_score_DSSP':seq_dist_SS}
temp_df = pd.DataFrame(temp_dict)
#print("temp_df........len",len(temp_df))

plt.figure(figsize=(20,8))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
sns.violinplot(data=temp_df, x='tm_score_min', y='SS_score_DSSP')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_SS_scale_DSSP_0_1_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_SS_scale_DSSP_0_1_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()






###### latters frequency

# H_alph = ['G', 'H', 'I']
# S_alph = ['E', 'B']
# L_alph = ['S', 'T', 'N']

df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_DSSP_seq_with_N_after_ck1-2-3.csv')

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
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/letter_distribution_SS_DSSP_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/letter_distribution_SS_DSSP_after_ck1-2-3.svg', format="svg")
plt.close()

