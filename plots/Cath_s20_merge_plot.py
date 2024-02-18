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



aa_df = pd.read_csv('/projects/globalscratch/alex/scope-nf/results4ballal/result.txt.gz',header=None)
aa_df.columns = ['name','Probab','E-value','Score','Aligned_cols','Identities','Similarity','Sum_probs','Template_Neff']
aa_df = aa_df[['name','E-value','Identities','Similarity']]


aa_df[['id','id2']]=aa_df['name'].str.split('-', expand=True)

##aa_df.columns = ['name','E-value','Identities','Similarity','id','id2']
##aa_df = aa_df[['Probab','E-value','Score','Aligned_cols','Identities','Similarity','Sum_probs','id','id2']]
aa_df['key_id']=['-'.join(sorted(combine)) for combine in zip(aa_df['id'], aa_df['id2'])]
aa_df['Identities'] = aa_df['Identities'].str.replace(r'\%', '')


tmp=aa_df.groupby('key_id')['Similarity']
##aa_df.assign(min_s=tmp.transform(min), max_s=tmp.transform(max))
aa_df['min_s'] = tmp.transform('min')
aa_df['max_s'] = tmp.transform('max')
tmp=aa_df.groupby('key_id')['Identities']
##aa_df.assign(min_s=tmp.transform(min), max_s=tmp.transform(max))
aa_df['min_I'] = tmp.transform('min')
aa_df['max_I'] = tmp.transform('max')

aa_df = aa_df.drop_duplicates(subset=['key_id'], keep='last')



aa_df = aa_df[['id', 'id2', 'key_id', 'min_s', 'max_s', 'min_I', 'max_I']]

aa_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/hh_result_s20.csv')





df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s20_results.csv')


df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-domain-description-v4_3_0.csv')
df_family =df_family[['id','Homol']]

df.columns = ['Unnamed: 0', 'id', 'domain_p1', 'ss_score', 'aa_score', 'p1_len',
       'p2_len', 'custom_local_idnt_p1_p2', 'custom_local_idnt_p2_p1',
       'key_id', 'TM_min', 'TM_max', '%id', 'E_value', 'qcov']


df = pd.merge(df, df_family ,how='left', on="id")
df.columns = ['Unnamed: 0', 'domain_p2', 'id', 'ss_score', 'aa_score', 'p1_len',
       'p2_len', 'custom_local_idnt_p1_p2', 'custom_local_idnt_p2_p1',
       'key_id', 'TM_min', 'TM_max', '%id', 'E_value', 'qcov', 'Homol_p2']

df = pd.merge(df, df_family ,how='left', on="id")
#print(df_merged.isnull().any())
df.columns = ['Unnamed: 0', 'domain_p2', 'domain_p1', 'ss_score', 'aa_score', 'p1_len',
       'p2_len', 'custom_local_idnt_p1_p2', 'custom_local_idnt_p2_p1',
       'key_id', 'TM_min', 'TM_max', '%id', 'E_value', 'qcov', 'Homol_p2',
       'Homol_p1']


df['cath_superFamily'] = np.where(df["Homol_p2"] == df["Homol_p1"], 1, 0)


df = df[['domain_p2', 'domain_p1', 'ss_score', 'aa_score',
       'p1_len', 'p2_len', 'custom_local_idnt_p1_p2',
       'custom_local_idnt_p2_p1', 'key_id', 'TM_min', 'TM_max', '%id',
       'E_value', 'qcov','cath_superFamily']]


# df_merged=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/cath_non_redundant_pdbs_s20_sequence_ForkPoolWorker-1.csv')
# df_merged = df_merged[~df_merged.AA_seq.isna()]
# df2 = df_merged.query('AA_seq.str.contains("\?")')
# df2_ids = df2['id'].tolist()
# df_merged = df_merged[~df_merged['id'].isin(df2_ids)]
#df_merged =df_merged.merge(df_family, left_on='id', right_on='id', how='inner')
# len(df_merged['Homol'].unique())
# df2_ids = df_merged['id'].tolist()
# len(df2_ids)
# df = df[df['domain_p2'].isin(df2_ids)]
# df = df[df['domain_p1'].isin(df2_ids)]
##df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s20_results.csv', index=False)

# df['E_value'] = df['E_value'].fillna(1)
# df['E_value'] = df['E_value'] + 0.000000001

# df['E_value'].isnull().values.any()
# df['E_value'] = np.log10(df['E_value'])


## x=df[df.E_value<0.05]


# df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s20_results.csv',usecols=['cath_superFamily','TM_min','ss_score'])
# df_dssp = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_s20_dssp_LEV_results.csv')

# df = df.merge(df_dssp,left_on='key_id',right_on='key_id')







df = pd.read_csv('/group/bioinf_protstr/All_alphaFoldSS_Jun2023/codes/v2/v3/v4/v5/SS/cath_s20/s20_merged.csv',usecols=['TM_min','LocalCust2_score','ss_score','cath_superFamily'])




## superfamily count
#count = (df['cath_superFamily'] == 0).sum()

########## AUC

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['TM_min'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...TM_min.",roc_auc)  ### 0.89816031807

plt.plot(fpr,tpr,color="blue",label="TM")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['ss_score'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...ss.",roc_auc,type(roc_auc))## 0.841362134
#create ROC curve
plt.plot(fpr,tpr,color="green",label="SS")

plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/tm_min_lev_roc_cath_s20.png', format="png")
plt.savefig('/group/bioinf/Ballal/plots/tm_min_lev_ss_roc_s20.svg', format="svg")
plt.close()





####### balanced_accuracy 

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['TM_min'], pos_label=1)

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
plt.savefig('/group/bioinf/Ballal/plots/TM_min_balance_accuracy_s20', format="svg")

plt.close()




## same / different superfamily plot

same_homo = df.query('cath_superFamily == 1') ['ss_score']
diff_homo = df.query('cath_superFamily == 0') ['ss_score']


same_homo_df = pd.DataFrame({'ss_score':same_homo.values, 'Catagory':'Same sfamily'})
diff_homo_df = pd.DataFrame({'ss_score':diff_homo.values, 'Catagory':'Different sfamily'})

frames = [same_homo_df, diff_homo_df]

same_homo_not_same_homo_df = pd.concat(frames)



g2 = sns.violinplot(data=same_homo_not_same_homo_df, y="ss_score", x="Catagory")
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf/Ballal/plots/SS_score_lev_same_different_superfamily_s20.svg', format="svg")
plt.close()









##TM_min

same_homo = df.query('cath_superFamily == 1') ['TM_min']
diff_homo = df.query('cath_superFamily == 0') ['TM_min']


same_homo_df = pd.DataFrame({'TM_min':same_homo.values, 'Catagory':'Same family'})
diff_homo_df = pd.DataFrame({'TM_min':diff_homo.values, 'Catagory':'Different family'})
 
frames = [same_homo_df, diff_homo_df]
same_homo_not_same_homo_df = pd.concat(frames)

g2 = sns.violinplot(data=same_homo_not_same_homo_df, y="TM_min", x="Catagory")
sns.despine()
g2.set(xlabel=None)
plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf/Ballal/plots/TM_min_same_different_superfamily_s20.svg', format="svg")
plt.close()


exit()

##

import pandas as pd


df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/cath_non_redundant_pdbs_s20_sequence_ForkPoolWorker-1.csv')

## id, SS_seq

df['length']  = df['SS_seq'].str.len()
print(df['length'].max())
print(df['length'].min())

sns.displot(df, x="length", kde=True)
plt.xlabel('Number of Residues')
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_all_cath_s20.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_all_cath_s20.svg', format="svg")
plt.close()





####


df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s20_results.csv',usecols=['cath_superFamily','key_id'])
hh_df =pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/hh_result_s20.csv',usecols=['min_I','key_id']) 

df = df.merge(hh_df,left_on='key_id',right_on='key_id')

fpr, tpr, thresholds = roc_curve(df["cath_superFamily"],df['min_I'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...min_I.",roc_auc)