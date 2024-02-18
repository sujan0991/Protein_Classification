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




df_lev = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_LEV_s40_results.csv')
df_dssp = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_s40_dssp_LEV_results.csv')


length_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/cath_non_redundant_pdbs_s40_sequence_ForkPoolWorker-1.csv')
length_df=length_df[~length_df.SS_seq.isna()]
length_df=length_df[~length_df.AA_seq.isna()]
length_df['ss_length'] = length_df['SS_seq'].str.len().astype(int)
length_df['aa_length'] = length_df['AA_seq'].str.len().astype(int)
length_df = length_df.query('ss_length == aa_length')
length_df = length_df[['id','ss_length']]

print(len(length_df))

df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv')
length_df =length_df.merge(df_family, left_on='id', right_on='id', how='inner')
length_df = length_df[['id','Homol']]

# len(length_df['Homol'].unique())

df_dssp['max_seq_len'] = df_dssp[["query_len", "subject_len"]].max(axis=1)
df_dssp['levIdent_SS'] = pd.to_numeric(df_dssp['SS_distance'])/ df_dssp['max_seq_len']
df_dssp['ss_score_dssp'] = abs(1-df_dssp["levIdent_SS"]).astype(np.float16)
df_dssp = df_dssp[['id', 'query_P','max_seq_len','ss_score_dssp']]

df_dssp = pd.merge(df_dssp, length_df ,how='left', on="id")
df_dssp=df_dssp[~df_dssp.Homol.isna()]
df_dssp.columns = ['subject', 'id', 'max_seq_len', 'ss_score_dssp', 'Homol_s']
df_dssp = pd.merge(df_dssp, length_df ,how='left', on="id")
df_dssp=df_dssp[~df_dssp.Homol.isna()]
df_dssp['cath_superFamily'] = np.where(df_dssp["Homol_s"] == df_dssp["Homol"], 1, 0)
df_dssp['key_id']=['-'.join(sorted(combine)) for combine in zip(df_dssp['id'], df_dssp['subject'])]
df_dssp = df_dssp[['max_seq_len', 'ss_score_dssp','cath_superFamily', 'key_id']]

##count = (df_dssp['cath_superFamily'] == 1).sum()

df_lev['key_id']=['-'.join(sorted(combine)) for combine in zip(df_lev['id'], df_lev['query_P'])]
df_lev = df_lev.merge(df_dssp,left_on='key_id',right_on='key_id')
df_dssp = pd.DataFrame()
df_lev['levIdent_SS'] = pd.to_numeric(df_lev['SS_distance'])/ df_lev['max_seq_len']
df_lev['ss_score'] = abs(1-df_lev["levIdent_SS"]).astype(np.float16)
df_lev['levIdent_AA'] = pd.to_numeric(df_lev['AA_distance'])/ df_lev['max_seq_len']
df_lev['aa_score'] = abs(1-df_lev["levIdent_AA"]).astype(np.float16)
df_lev = df_lev[['key_id','ss_score_dssp', 'cath_superFamily','ss_score', 'aa_score']]





df_tm = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_s40_results.csv')

df_tm.columns = ['Unnamed: 0', 'subject', 'query', 'TM_query', 'TM_subject']
df_tm['key_id']=['-'.join(sorted(combine)) for combine in zip(df_tm['subject'], df_tm['query'])]

df_lev = df_lev.merge(df_tm,left_on='key_id',right_on='key_id')

df_tm = pd.DataFrame()

df_lev['TM_max'] = df_lev[["TM_query", "TM_subject"]].max(axis=1)
df_lev['TM_min'] = df_lev[["TM_query", "TM_subject"]].min(axis=1)



df_lev = df_lev[['subject', 'query','key_id', 'ss_score_dssp', 'cath_superFamily', 'ss_score', 'aa_score','TM_max','TM_min']]






df_blast = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_non_redundant_pdbs_s40_blast.csv')
df_blast.columns = ['query','subject','%id','alignment_length','mismatches','gap_openings','query_start','query_end','subject_start','subject_end','E_value','bit_score','qcov']
df_blast['key_id']=['-'.join(sorted(combine)) for combine in zip(df_blast['subject'], df_blast['query'])]
df_blast.drop_duplicates(subset='key_id', keep="last", inplace=True)
df_blast = df_blast[['%id','E_value','qcov', 'key_id']]
df_lev = df_lev.merge(df_blast,on='key_id',how='left')


##df_lev.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s40_results.csv', index=False)




#df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s40_results.csv')


# Use `.loc` to modify the DataFrame
df_lev.loc[(df_lev['E_value'] > 100) | (df_lev['E_value'].isna()), 'E_value'] = 100






########## AUC

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['TM_min'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...TM_min.",roc_auc)  ### 0.9307323

plt.plot(fpr,tpr,color="blue",label="TM")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')


fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['ss_score'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))## 0.7560


#create ROC curve
plt.plot(fpr,tpr,color="green",label="SS")
#plt.plot([], [], ' ', label= round(roc_auc, 2))

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['aa_score'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc....aa",roc_auc)## 0.5864395

plt.plot(fpr,tpr,color="orange",label="AA")



plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/tm_min_lev_roc_cath_s20.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s40_plots/tm_min_lev_roc_cath_s40.svg', format="svg")
plt.close()





fpr, tpr, thresholds = roc_curve(df_lev["cath_superFamily"],df_lev['E_value'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...E_value.",roc_auc) ### 0.48302082775

fpr, tpr, thresholds = roc_curve(df_lev["cath_superFamily"],df_lev['ss_score_dssp'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...ss_score_dssp.",roc_auc) ### 0.722530019


####### balanced_accuracy 

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['aa_score'], pos_label=1)

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


plt.figtext(0.41, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12) ## for aa x=0.30,ss=0.55
plt.figtext(0.08, 0.73, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)## for aa y=0.59,ss= 0.78

#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/TM_min_balance_accuracy.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s40_plots/aa_score_balance_accuracy.svg', format="svg")

plt.close()


