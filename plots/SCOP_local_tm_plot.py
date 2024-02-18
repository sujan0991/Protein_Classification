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



lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_HH_blast_results.csv',
                        usecols=['is_same_scope_Superfamily','SS_score','key_id','tm_min'])
local_df = pd.read_csv('/group/bioinf/Ballal/scop_local_custom_weight_results_full.csv')
local_df = local_df.query("`query_id` != 'query_id'")
local_df.to_csv('/group/bioinf/Ballal/scop_local_custom_weight_results_full.csv')

lev_tm_df.columns = ['is_same_scope_Superfamily','SS_score_lev','key_id','tm_min']
local_df.columns = ['query_id', 'id', 'local_sim', 'local_alignLength', 'length_query','length_target']
local_df['key_id']=['-'.join(sorted(combine)) for combine in zip(local_df['query_id'], local_df['id'])]
lev_tm_df['max_seq_len'] = lev_tm_df[["length_query", "length_target"]].max(axis=1)
lev_tm_df['local_sim'] = lev_tm_df['local_sim'].astype(float)
lev_tm_df['LocalCust2_score'] = lev_tm_df['local_sim']/ lev_tm_df['max_seq_len']



## local_df.isnull().values.any()

lev_tm_df = lev_tm_df.merge(local_df,left_on='key_id',right_on='key_id')



lev_tm_df.to_csv('/group/bioinf/Ballal/scop_lev_ss_tm_local_custom_weight_results_full.csv')

lev_tm_df = pd.read_csv('/group/bioinf/Ballal/scop_lev_ss_tm_local_custom_weight_results_full.csv')

## AUC


fpr, tpr, thresholds = roc_curve(lev_tm_df["is_same_scope_Superfamily"],lev_tm_df['tm_min'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...TM_min.",roc_auc)  ### 0.9554
plt.plot(fpr,tpr,color="blue",label="TM")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')

fpr, tpr, thresholds = roc_curve(lev_tm_df["is_same_scope_Superfamily"],lev_tm_df['SS_score_lev'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))## 0.89727
plt.plot(fpr,tpr,color="green",label="SS")
plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')
plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/tm_min_lev_roc_cath_s20.png', format="png")
plt.savefig('/group/bioinf/Ballal/plots/tm_min_lev_ss_roc_scop.svg', format="svg")
plt.close()



####### balanced_accuracy 

fpr, tpr, thresholds = roc_curve(lev_tm_df["is_same_scope_Superfamily"],lev_tm_df['LocalCust2_score'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2


best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

plt.axis([0, 1, 0, 1]) ## comment it for HHBlits
plt.plot(thresholds,tpr, color="blue",label="tpr" )
#plt.plot(thresholds,fpr, color="purple",label="fpr" )
plt.plot(thresholds,tnr, color="green", label="tnr")
#plt.plot(thresholds,fnr, color="orange",label="fnr" )
plt.plot(thresholds,balanced_accuracy, color="black",label="ba")
plt.legend(loc="lower right")


plt.figtext(0.41, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12) ## for aa x=0.30,ss=0.55
plt.figtext(0.08, 0.73, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)## for aa y=0.59,ss= 0.78

#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/TM_min_balance_accuracy.png', format="png")
plt.savefig('/group/bioinf/Ballal/plots/ss_local_balance_accuracy_scop', format="svg")

plt.close()




## same / different superfamily plot

same_homo = lev_tm_df.query('is_same_scope_Superfamily == 1') ['SS_score_lev']
diff_homo = lev_tm_df.query('is_same_scope_Superfamily == 0') ['SS_score_lev']


same_homo_df = pd.DataFrame({'ss_score':same_homo.values, 'Catagory':'Same family'})
diff_homo_df = pd.DataFrame({'ss_score':diff_homo.values, 'Catagory':'Different family'})

frames = [same_homo_df, diff_homo_df]

same_homo_not_same_homo_df = pd.concat(frames)



g2 = sns.violinplot(data=same_homo_not_same_homo_df, y="ss_score", x="Catagory")
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf/Ballal/plots/SS_score_lev_same_different_superfamily_scop.svg', format="svg")
plt.close()









##TM_min

same_homo = lev_tm_df.query('is_same_scope_Superfamily == 1') ['tm_min']
diff_homo = lev_tm_df.query('is_same_scope_Superfamily == 0') ['tm_min']


same_homo_df = pd.DataFrame({'TM_min':same_homo.values, 'Catagory':'Same family'})
diff_homo_df = pd.DataFrame({'TM_min':diff_homo.values, 'Catagory':'Different family'})
 
frames = [same_homo_df, diff_homo_df]
same_homo_not_same_homo_df = pd.concat(frames)

g2 = sns.violinplot(data=same_homo_not_same_homo_df, y="TM_min", x="Catagory")
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf/Ballal/plots/TM_min_same_different_superfamily_scop.svg', format="svg")
plt.close()

