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


df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_results.csv')
df = df[df.AA_distance != "AA_distance"]


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
family_df["scope_superFamily"] = family_df['scop'].str.rpartition('.')[0]

len(family_df['scope_superFamily'].unique())

length_df =length_df.merge(family_df, left_on='id', right_on='id', how='inner')

###########
##length_df =length_df[['id','ss_length','scop']]

df = pd.merge(df, length_df ,how='left', on="id")

##df = df[['id', 'query_P', 'AA_distance', 'SS_distance', 'ss_length', 'scop']]
df=df[~df.scop.isna()]
df.columns = ['protein1', 'id', 'AA_distance', 'SS_distance', 'protein1_length', 'scop1']
df = pd.merge(df, length_df ,how='left', on="id")
df=df[~df.scop.isna()]
df.columns = ['protein1','protein2','AA_distance', 'SS_distance', 'protein1_length', 'scop1','protein2_length', 'scop2']


df['max_seq_len'] = df[["protein1_length", "protein2_length"]].max(axis=1)
df['levIdent_SS'] = pd.to_numeric(df['SS_distance'])/ df['max_seq_len']
df['SS_score'] = abs(1-df["levIdent_SS"]).astype(np.float16)
df['levIdent_AA'] = pd.to_numeric(df['AA_distance'])/ df['max_seq_len']
df['AA_score'] = abs(1-df["levIdent_AA"]).astype(np.float16)
df['scope_Family'] = np.where(df["scop1"] == df["scop2"], 1, 0)


df["scope_superFamily1"] = df['scop1'].str.rpartition('.')[0]
df["scope_superFamily2"] = df['scop2'].str.rpartition('.')[0]

df['is_same_scope_Superfamily'] = np.where(df["scope_superFamily1"] == df["scope_superFamily2"], 1, 0)


# df['class1'] = df['scop1'].astype(str).str[0]
# df['class2'] = df['scop2'].astype(str).str[0]



df = df[['protein1','protein2', 'scope_Family','is_same_scope_Superfamily','SS_score','AA_score']]

df['key_id']=['-'.join(sorted(combine)) for combine in zip(df['protein1'], df['protein2'])]


tm_df = pd.read_csv('/group/bioinf_protstr/article_briefings/data/scope/structure/structure.csv')
tm_df["tm_min"] = tm_df[["tm_query", "tm_subject"]].min(axis=1)
tm_df["tm_max"] = tm_df[["tm_query", "tm_subject"]].max(axis=1)
tm_df = tm_df[['query', 'subject', 'tm_min', 'tm_max']]
tm_df['key_id']=['-'.join(sorted(combine)) for combine in zip(tm_df['query'], tm_df['subject'])]


df = df.merge(tm_df,left_on='key_id',right_on='key_id')

df = df[['protein1', 'protein2', 'scope_Family','is_same_scope_Superfamily', 'SS_score', 'AA_score', 'key_id','tm_min', 'tm_max']]


aa_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/hhalign/parsed_hhalign_results.csv',header=None)
tmp=aa_df.groupby('key_id')['Identities']
##aa_df.assign(min_s=tmp.transform(min), max_s=tmp.transform(max))
aa_df['min_I'] = tmp.transform('min')
aa_df['max_I'] = tmp.transform('max')
aa_df = aa_df[['id', 'id2', 'key_id','min_I', 'max_I']]
aa_df = aa_df.drop_duplicates(subset=['key_id'], keep='last')
df = df.merge(aa_df,left_on='key_id',right_on='key_id')


blast_df = pd.read_csv('/group/bioinf_protstr/article_briefings/data/scope/sequence/blast.csv')
blast_df['blast_identity'] = (blast_df['%id']/100)* blast_df['qcov']
blast_df = blast_df[['query','subject','blast_identity','qcov','E_value']]
blast_df.columns = ['id','subject','blast_identity','qcov','E_value']

blast_df['key_id']=['-'.join(sorted(combine)) for combine in zip(blast_df['id'], blast_df['subject'])]
tmp=blast_df.groupby('key_id')['blast_identity']
blast_df['blast_identity_min'] = tmp.transform('min')
blast_df['blast_identity_max'] = tmp.transform('max')

tmp=blast_df.groupby('key_id')['qcov']
blast_df['qcov_min'] = tmp.transform('min')
blast_df['qcov_max'] = tmp.transform('max')

tmp=blast_df.groupby('key_id')['E_value']
blast_df['E_value_min'] = tmp.transform('min')
blast_df['E_value_max'] = tmp.transform('max')


blast_df.drop_duplicates(subset='key_id', keep="last", inplace=True)

df = pd.merge(df, blast_df ,how='left', on="key_id")

##df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_HH_blast_results.csv')

# df= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_HH_blast_results.csv')


##df.isnull().values.any()


#############

### for numbers

## min / same / diff

same = df.query('is_same_scope_Superfamily == 1 & SS_score <= 0.5')
print('>=0.5',len(same))


same = df.query('is_same_scope_Superfamily == 1 & SS_score > 0.5 & SS_score < 0.75')
print('2-3',len(same))


same = df.query('is_same_scope_Superfamily == 1 & SS_score >= 0.75')
print('< 0.5',len(same))


same = df.query('is_same_scope_Superfamily == 1 & tm_max < 0.5')
print('< 0.5',len(same))
same = df.query('is_same_scope_Superfamily == 1 & tm_max >= 0.5')
print('>= 0.5',len(same))


## max

same = df.query('tm_max >= 0.5 & SS_score <= 0.5')
print('3',len(same))






### plots

###ss

same_homo = df.query('scope_Family == 1') ['SS_score']
diff_homo = df.query('scope_Family == 0') ['SS_score']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.7:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same family:total',len(same_homo), 'above point 0.7', above_point_5 ,'below point 0.7', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.7:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different family:total',len(diff_homo), ' above point 0.7', above_point_5 ,'below point 0.7', below_point_5)
 
same_homo_df = pd.DataFrame({'SS_score':same_homo.values, 'Catagory':'Same family'})
diff_homo_df = pd.DataFrame({'SS_score':diff_homo.values, 'Catagory':'Different family'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

same_homo_not_same_homo_df['SS_score'].isnull().values.any()

g2 = sns.violinplot(data=same_homo_not_same_homo_df, y="SS_score", x="Catagory")
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_family_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_family_scope.svg', format="svg")
plt.close()


## is_same_scope_Superfamily

same_homo = df.query('is_same_scope_Superfamily == 1') ['SS_score']
diff_homo = df.query('is_same_scope_Superfamily == 0') ['SS_score']

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
 
same_homo_df = pd.DataFrame({'SS_score':same_homo.values, 'Catagory':'Same superfamily'})
diff_homo_df = pd.DataFrame({'SS_score':diff_homo.values, 'Catagory':'Different superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

same_homo_not_same_homo_df['SS_score'].isnull().values.any()

g2 = sns.violinplot(data=same_homo_not_same_homo_df, y="SS_score", x="Catagory")
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_superfamily_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/SS_score_superfamily_scope.svg', format="svg")
plt.close()



## AA

same_homo = df.query('is_same_scope_Superfamily == 1') ['AA_score']
diff_homo = df.query('is_same_scope_Superfamily == 0') ['AA_score']

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

g2 = sns.violinplot( y=same_homo_not_same_homo_df["AA_score"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/AA_score_superfamily_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/AA_score_superfamily_scope_.svg', format="svg")
plt.close()



#######

same_homo = df.query('is_same_scope_Superfamily == 1') ['tm_max']
diff_homo = df.query('is_same_scope_Superfamily == 0') ['tm_max']

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
 
same_homo_df = pd.DataFrame({'tm_min':same_homo.values, 'Catagory':'Same superfamily'})
diff_homo_df = pd.DataFrame({'tm_min':diff_homo.values, 'Catagory':'Different superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

same_homo_not_same_homo_df['tm_min'].isnull().values.any()

g2 = sns.violinplot(data=same_homo_not_same_homo_df, y="tm_min", x="Catagory")
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_score_superfamily_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/TM_max_score_superfamily_scope.svg', format="svg")
plt.close()





df2=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_LEV_TM_results.csv',usecols=['is_same_scope_Superfamily','tm_min','SS_score_DSSP'])


###### auc

fpr, tpr, thresholds = roc_curve(df2["is_same_scope_Superfamily"],df2['tm_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...tm.min",roc_auc)

plt.plot(fpr,tpr,color="blue",label="TM_min")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')

fpr, tpr, thresholds = roc_curve(df2["is_same_scope_Superfamily"],df2['SS_score_DSSP'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...SS_score_DSSP.",roc_auc,type(roc_auc))


#create ROC curve
plt.plot(fpr,tpr,color="green",label="DSSP")


#plt.plot([], [], ' ', label= round(roc_auc, 2))


# fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['min_I'], pos_label=1)
# print("fpr, tpr len",len(fpr),len(tpr))
# roc_auc = auc(fpr, tpr)
# print("roc_auc....min_I",roc_auc)

# plt.plot(fpr,tpr,color="orange",label="min_I")


plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_min_lev_min_I_roc_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_min_dssp_roc_scope.svg', format="svg")
plt.close()


###

df['min_I'] = df['min_I'].astype(float)

fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['min_I'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc....min_I",roc_auc)

df['max_I'] = df['max_I'].astype(float)

fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['max_I'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc....max_I",roc_auc)



df['blast_identity_min'] = df['blast_identity_min'].fillna(0)
df['blast_identity_min'].isnull().values.any()
df['blast_identity_max'] = df['blast_identity_max'].fillna(0)
df['blast_identity_max'].isnull().values.any()

df['E_value_min'] = df['E_value_min'].fillna(1)
df['E_value_max'] = df['E_value_max'].fillna(1)
df['E_value_min'] = df['E_value_min'] + 0.000000001
df['E_value_max'] = df['E_value_max'] + 0.000000001

df['E_value_min'].isnull().values.any()
df['E_value_max'].isnull().values.any()

df['E_value_min'] = np.log10(df['E_value_min'])
df['E_value_max'] = np.log10(df['E_value_max'])



fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['tm_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc....tm_min",roc_auc)


#### balance accuracy


fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['min_I'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2

best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
#balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

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

# plt.axvline(x=thresholds[(thresholds>=0.5)&(thresholds<=0.5)],ymax=round(balanced_accuracy_05[0],2), color="darkgrey",linestyle='--')
# plt.axhline(y=round(balanced_accuracy_05[0],2),xmax=thresholds[(thresholds>=0.5)&(thresholds<=0.5)], color="darkgrey",linestyle='--')
# plt.figtext(0.08, 0.76, round(balanced_accuracy_05[0],2), wrap=True, horizontalalignment='center', fontsize=12) #### ss = 0.74



#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/max_I_superfamily_balanced_accuracy_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/min_I_superfamily_balanced_accuracy_scope.svg', format="svg")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/aa_superfamily_balanced_accuracy_scope.svg', format="svg")

plt.close()






##

seq_dist_SS = df['SS_score'].to_list()
tm_score_min = df["tm_min"].to_list()

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
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_SS_scale_0_1_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_SS_scale_0_1_scope.svg', format="svg")
plt.close()






######### distribution


length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_ForkPoolWorker-1.csv')
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')

length_df['class'] = length_df['scop'].astype(str).str[0]


## %NUMBER OF RESI

sns.displot(length_df, x="aa_length", kde=True)
plt.xlabel('Number of Residues')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_scope.svg', format="svg")
plt.close()



### class distribution


cl_count = length_df['class'].value_counts()

df = pd.DataFrame({'class':cl_count.index, 'count':cl_count.values})


plt.figure(figsize=(10,8))
sns.barplot(data=df, x="class", y="count")
plt.ylabel('Count')
plt.xlabel('')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/class_distribution_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/class_distribution_scope.svg', format="svg")
plt.close()



### %PYMOL _Letter_frequency

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
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/letter_distribution_SS_scope.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/letter_distribution_SS_scope.svg', format="svg")
plt.close()



