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


df = pd.read_csv('/group/bioinf_protstr/Ballal/Foldseek_results/foldseekTM.txt',header=None)

df.columns = ['all']

df = pd.concat([df['all'], df['all'].str.split('\t', expand = True)], axis = 1)

df = df[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]]

df[0] = df[0].str.replace(r'.pdb', '')
df[1] = df[1].str.replace(r'.pdb', '')
df['key_id']=['-'.join(sorted(combine)) for combine in zip(df[0], df[1])]

df = df[[2, 'key_id']]

df.columns = ['foldseek_identity', 'key_id']

tmp=df.groupby('key_id')['foldseek_identity']
df['min_F_I'] = tmp.transform('min')
df['max_F_I'] = tmp.transform('max')



df.drop_duplicates(subset='key_id', keep="last", inplace=True)



ss_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_LEV_TM_results.csv')

##df = df.merge(ss_df,left_on='key_id',right_on='key_id')
df =df.merge(ss_df, left_on='key_id', right_on='key_id', how='inner')

hh_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_HH_blast_results.csv')
hh_df = hh_df[['key_id','min_I', 'max_I','blast_identity_min', 'blast_identity_max']]
df =df.merge(hh_df, left_on='key_id', right_on='key_id', how='inner')


df = df [['foldseek_identity', 'key_id', 'min_F_I', 'max_F_I', 'protein1_x', 'protein2_x', 
          'is_same_scope_Superfamily_x', 'SS_score_DSSP_x', 'tm_min_x', 'tm_max_x', 'AA_score_x', 'SS_score_x',
          'min_I', 'max_I','blast_identity_min', 'blast_identity_max']]

df.columns = ['foldseek_identity', 'key_id', 'min_F_I', 'max_F_I', 'protein1', 'protein2', 'is_same_scope_Superfamily', 
              'SS_score_DSSP', 'tm_min', 'tm_max', 'AA_score', 'SS_score', 'min_I', 'max_I', 'blast_identity_min', 'blast_identity_max']

## auc

df['min_F_I'] = df['min_F_I'].astype(float)
df['max_F_I'] = df['max_F_I'].astype(float)
df['foldseek_identity'] = df['foldseek_identity'].astype(float)
##df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/foldseek_merged_result.csv')



fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['max_F_I'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...max_F_I.",roc_auc,type(roc_auc))


fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['SS_score'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))

fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['AA_score'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...AA_score.",roc_auc,type(roc_auc))


fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['tm_max'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...tm.max",roc_auc)

fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['SS_score_DSSP'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...SS_score_DSSP",roc_auc)

fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['min_I'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...min_I.",roc_auc,type(roc_auc))

fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['max_I'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...max_I.",roc_auc,type(roc_auc))


df['blast_identity_min'] = df['blast_identity_min'].fillna(0)
df['blast_identity_min'].isnull().values.any()

df['blast_identity_max'] = df['blast_identity_max'].fillna(0)
df['blast_identity_max'].isnull().values.any()


fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['blast_identity_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity_min.",roc_auc)


fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['blast_identity_max'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity_max.",roc_auc)



### balance accuracy
##'min_F_I', 'max_F_I''min_I', 'max_I', 'blast_identity_min', 'blast_identity_max'
fpr, tpr, thresholds = roc_curve(df["is_same_scope_Superfamily"],df['blast_identity_max'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2


best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

#plt.axis([0, 1, 0, 1])
plt.plot(thresholds,tpr, color="blue",label="tpr" )
plt.plot(thresholds,fpr, color="purple",label="fpr" )
plt.plot(thresholds,tnr, color="green", label="tnr")
plt.plot(thresholds,fnr, color="orange",label="fnr" )
plt.plot(thresholds,balanced_accuracy, color="black",label="ba")
plt.legend(loc="lower right")

plt.figtext(0.41, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12) ## for aa x=0.30,ss=0.55
plt.figtext(0.08, 0.73, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)## for aa y=0.59,ss= 0.78

plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_max_superfamily_balanced_accuracy_scope_foldseek.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/blast_identity_max_superfamily_balanced_accuracy_scope_foldseek.svg', format="svg")

plt.close()

