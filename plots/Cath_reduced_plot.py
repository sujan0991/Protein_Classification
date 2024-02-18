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


df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_LEV_reduced_results.csv')

length_df = pd.read_csv('/group/bioinf_protstr/toAli/cath_SS_seq_reduced.csv', usecols=['id','SS_seq_reduced'])
length_df=length_df[~length_df.SS_seq_reduced.isna()]
length_df['ss_length'] = length_df['SS_seq_reduced'].str.len().astype(int)

df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')

df = pd.merge(df, length_df ,how='left', on="id")
df = df[['id', 'query_P','SS_distance', 'ss_length', 'Homol']]
df.columns = ['domain1','id', 'SS_distance', 'domain1_length','Homol1']
df = pd.merge(df, length_df ,how='left', on="id")
df = df[['domain1', 'id', 'SS_distance', 'domain1_length', 'Homol1','ss_length','Homol']]
df.columns = ['domain1','domain2','SS_distance', 'domain1_length', 'Homol1','domain2_length', 'Homol2']


df['max_seq_len'] = df[["domain1_length", "domain2_length"]].max(axis=1)
df['levIdent_SS'] = pd.to_numeric(df['SS_distance'])/ df['max_seq_len']
df['SS_score_reduced'] = abs(1-df["levIdent_SS"]).astype(np.float16)
df['cath_superFamily'] = np.where(df["Homol1"] == df["Homol2"], 1, 0)
df = df[['domain1','domain2', 'cath_superFamily','SS_score_reduced']]

df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_ss_score_reduced_results.csv')


lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_lev_tm_HH_blast_blast_ss_results.csv',usecols=['cath_superFamily','TM_min','SS_score'])


## AUC

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['TM_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...tm.max",roc_auc)

plt.plot(fpr,tpr,color="blue",label="TM_max")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['SS_score_reduced'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))

#create ROC curve
plt.plot(fpr,tpr,color="green",label="SS")



plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/tm_max_lev_min_I_roc_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/thesis_presentation_plots/auc/svg/tm_min_lev_reduced_roc_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()









