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



blast_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/BLAST_CATH/blast.csv')
print(len(blast_df))  ## 3,877,669
blast_df['blast_identity'] = (blast_df['%id']/100)* blast_df['qcov']
blast_df = blast_df[['query','subject','blast_identity','qcov','E_value']]
blast_df.columns = ['id','subject','blast_identity','qcov','E_value']

df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-domain-description-v4_3_0.csv')
blast_df =blast_df.merge(df_family, left_on='id', right_on='id', how='inner')

##['id', 'subject', 'blast_identity', 'index', 'Class', 'Arch', 'Topol', 'Homol']
blast_df = blast_df[['id', 'subject', 'blast_identity','qcov','E_value', 'Class', 'Arch', 'Topol', 'Homol']]
blast_df.columns = ['domain1','id','blast_identity','qcov','E_value','Class1', 'Arch1', 'Topol1', 'Homol1']
blast_df =blast_df.merge(df_family, left_on='id', right_on='id', how='inner')
##['domain1', 'id', 'blast_identity', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'index', 'Class', 'Arch', 'Topol', 'Homol']
blast_df = blast_df[['domain1', 'id', 'blast_identity','qcov','E_value', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'Class', 'Arch', 'Topol', 'Homol']]
blast_df.columns = ['domain1', 'domain2', 'blast_identity','qcov','E_value', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'Class2', 'Arch2', 'Topol2', 'Homol2']

blast_df['cath_superFamily'] = np.where(blast_df["Homol1"] == blast_df["Homol2"], 1, 0)

blast_df = blast_df[['domain1', 'domain2', 'blast_identity','cath_superFamily','qcov','E_value']]


blast_df['blast_identity*qcov'] = blast_df['blast_identity']* blast_df['qcov']


#### blast identity

same_homo = blast_df.query('cath_superFamily == 1') ['blast_identity']
diff_homo = blast_df.query('cath_superFamily == 0') ['blast_identity']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 50.0:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 50', above_point_5 ,'below point 50', below_point_5)
##lev for same superfamily:total 3687494 above point 50 2027409 below point 50 1660085


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 50.0:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 0.5', above_point_5 ,'below point 0.5', below_point_5)
## lev for different superfamily:total 190175  above point 0.5 708 below point 0.5 189467

same_homo_df = pd.DataFrame({'blast_identity':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/png/blast_identity_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/svg/blast_identity_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


##lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_TM_results_all_safe_cath_after_ck1-2-3.csv')


##### blast_identity*qcov


same_homo = blast_df.query('cath_superFamily == 1') ['blast_identity*qcov']
diff_homo = blast_df.query('cath_superFamily == 0') ['blast_identity*qcov']


same_homo_df = pd.DataFrame({'blast_identity*qcov':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity*qcov':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity*qcov"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/png/blast_identity_qcov_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/svg/blast_identity_qcov_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()


##### blast_ e value


same_homo = blast_df.query('cath_superFamily == 1') ['E_value']
diff_homo = blast_df.query('cath_superFamily == 0') ['E_value']


same_homo_df = pd.DataFrame({'E_value':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'E_value':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["E_value"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/png/blast_E_value_super_safe_cath_after_ck1-2-3.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/svg/blast_E_value_super_safe_cath_after_ck1-2-3.svg', format="svg")
plt.close()
