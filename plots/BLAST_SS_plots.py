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


#### cat *txt > /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_ss_BLAST_result_all_replaced_S_H_L.txt

df=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_ss_BLAST_result_all_replaced_S_H_L.txt', header=None, sep='\t')
df.columns=['id', 'subject', 'blast_identity','qcov','qlen' ,'slen' ,'length','bitscore','e_value','seq1','seq2']
df['key_id']=['-'.join(sorted(combine)) for combine in zip(df['id'], df['subject'])]
df=df.sort_values(['key_id', 'e_value'], ascending=[True, True])
##df[['key_id', 'e_value']]

df.drop_duplicates(subset='key_id', keep="first", inplace=True)
df['blast_identity_qcov'] = (df['blast_identity']* df['qcov'])/100

df=df[['key_id','blast_identity', 'qcov', 'blast_identity_qcov','e_value']]

lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_TM_results_all_safe_cath_after_ck1-2-3.csv')

lev_tm_df= lev_tm_df.merge(df, on='key_id', how='left')
lev_tm_df['blast_identity_qcov'] = lev_tm_df['blast_identity_qcov'].fillna(0)
lev_tm_df['blast_identity'] = lev_tm_df['blast_identity'].fillna(0)
lev_tm_df['qcov'] = lev_tm_df['qcov'].fillna(0)
lev_tm_df['e_value'] = lev_tm_df['e_value'].fillna(1)
##lev_tm_df['e_value'].replace(0, 1, inplace=True)


lev_tm_df['logEval'] = np.log10((lev_tm_df.e_value)+0.00000000000000001)

######## plots


### AUC

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['blast_identity_qcov'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity_qcov.",roc_auc)
##roc_auc...blast_identity_qcov. 0.6006362580343847

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['blast_identity'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity.",roc_auc)
##roc_auc...blast_identity. 0.6011940596645293

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['qcov'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...qcov.",roc_auc)
##roc_auc...qcov. 0.5996705557721058

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['logEval'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...logEval.",roc_auc)
##roc_auc...logEval. 0.4354143887764534



## blast_identity_qcov
same_homo = lev_tm_df.query('cath_superFamily == 1') ['blast_identity_qcov']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['blast_identity_qcov']
same_homo_df = pd.DataFrame({'blast_identity_qcov':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity_qcov':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity_qcov"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/blast_identity_qcov_ss_same_diff_superF_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()


### blast_identity

same_homo = lev_tm_df.query('cath_superFamily == 1') ['blast_identity']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['blast_identity']
same_homo_df = pd.DataFrame({'blast_identity':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/blast_identity_ss_same_diff_superF_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()


## qcov

same_homo = lev_tm_df.query('cath_superFamily == 1') ['qcov']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['qcov']
same_homo_df = pd.DataFrame({'qcov':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'qcov':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["qcov"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/qcov_ss_same_diff_superF_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()


## logEval

same_homo = lev_tm_df.query('cath_superFamily == 1') ['logEval']
diff_homo = lev_tm_df.query('cath_superFamily == 0') ['logEval']
same_homo_df = pd.DataFrame({'logEval':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'logEval':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["logEval"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/logEval_ss_same_diff_superF_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()





#### balanced_accuracy

fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['logEval'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2

best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]

plt.axis([0, 1, 0, 1])
plt.plot(thresholds,tpr, color="blue",label="tpr" )
plt.plot(thresholds,fpr, color="purple",label="fpr" )
plt.plot(thresholds,tnr, color="green", label="tnr")
plt.plot(thresholds,fnr, color="orange",label="fnr" )
plt.plot(thresholds,balanced_accuracy, color="black",label="ba")
plt.legend(loc="lower right")

plt.figtext(0.41, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12) #
plt.figtext(0.08, 0.73, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/logEval_ss_balance_ac.png', format="png")
plt.close()


#### nonzero values

lev_tm_df_NonZero = lev_tm_df[lev_tm_df.blast_identity_qcov!=0]
#lev_tm_df_NonZero['logEval'] = np.log10((lev_tm_df_NonZero.e_value)+0.00000000000000001)
######## plots


### AUC

fpr, tpr, thresholds = roc_curve(lev_tm_df_NonZero["cath_superFamily"],lev_tm_df_NonZero['blast_identity_qcov'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity_qcov.",roc_auc)
##roc_auc...blast_identity_qcov. 0.7147559520485792

fpr, tpr, thresholds = roc_curve(lev_tm_df_NonZero["cath_superFamily"],lev_tm_df_NonZero['blast_identity'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...blast_identity.",roc_auc)
##roc_auc...blast_identity. 0.7901575153191764

fpr, tpr, thresholds = roc_curve(lev_tm_df_NonZero["cath_superFamily"],lev_tm_df_NonZero['qcov'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...qcov.",roc_auc)
##roc_auc...qcov. 0.5842021384058506

fpr, tpr, thresholds = roc_curve(lev_tm_df_NonZero["cath_superFamily"],lev_tm_df_NonZero['logEval'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc...logEval.",roc_auc)
##roc_auc...logEval. 0.38840107388251754


## 
fpr, tpr, thresholds = roc_curve(lev_tm_df_NonZero["cath_superFamily"],lev_tm_df_NonZero['qcov'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2

best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]

##plt.axis([0, 1, 0, 1])
plt.plot(thresholds,tpr, color="blue",label="tpr" )
plt.plot(thresholds,fpr, color="purple",label="fpr" )
plt.plot(thresholds,tnr, color="green", label="tnr")
plt.plot(thresholds,fnr, color="orange",label="fnr" )
plt.plot(thresholds,balanced_accuracy, color="black",label="ba")
plt.legend(loc="lower right")

plt.figtext(0.41, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12) #
plt.figtext(0.08, 0.73, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/qcov_ss_balance_ac_only_blast_data.png', format="png")
plt.close()


#########

same_homo = lev_tm_df_NonZero.query('cath_superFamily == 1') ['blast_identity_qcov']
diff_homo = lev_tm_df_NonZero.query('cath_superFamily == 0') ['blast_identity_qcov']
same_homo_df = pd.DataFrame({'blast_identity_qcov':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity_qcov':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity_qcov"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/blast_identity_qcov_ss_same_diff_superF_only_blast_data_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()


### blast_identity 
same_homo = lev_tm_df_NonZero.query('cath_superFamily == 1') ['blast_identity']
diff_homo = lev_tm_df_NonZero.query('cath_superFamily == 0') ['blast_identity']
same_homo_df = pd.DataFrame({'blast_identity':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'blast_identity':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/blast_identity_ss_same_diff_superF_only_blast_data_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()

### qcov 
same_homo = lev_tm_df_NonZero.query('cath_superFamily == 1') ['qcov']
diff_homo = lev_tm_df_NonZero.query('cath_superFamily == 0') ['qcov']
same_homo_df = pd.DataFrame({'qcov':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'qcov':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["qcov"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/qcov_ss_same_diff_superF_only_blast_data_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()



## logEval
same_homo = lev_tm_df_NonZero.query('cath_superFamily == 1') ['logEval']
diff_homo = lev_tm_df_NonZero.query('cath_superFamily == 0') ['logEval']
same_homo_df = pd.DataFrame({'logEval':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'logEval':diff_homo.values, 'Catagory':'Different Superfamily'})
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["logEval"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#plt.ylim(0, 1)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/logEval_ss_same_diff_superF_only_blast_data_v.png', format="png")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/y.svg', format="svg")
plt.close()



####################

# blast_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/BLAST_CATH/blast_SS.csv')
# blast_df['blast_identity_ss'] = (blast_df['%id']/100)* blast_df['qcov']
# blast_df = blast_df[['query','subject','blast_identity_ss','qcov','E_value']]
# blast_df.columns = ['id','subject','blast_identity_ss','qcov_ss','E_value_ss']

# df_family= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath-domain-description-v4_3_0.csv')
# blast_df =blast_df.merge(df_family, left_on='id', right_on='id', how='inner')

# blast_df = blast_df[['id', 'subject', 'blast_identity','qcov','E_value', 'Class', 'Arch', 'Topol', 'Homol']]
# blast_df.columns = ['domain1','id','blast_identity','qcov','E_value','Class1', 'Arch1', 'Topol1', 'Homol1']
# blast_df =blast_df.merge(df_family, left_on='id', right_on='id', how='inner')

# blast_df = blast_df[['domain1', 'id', 'blast_identity','qcov','E_value', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'Class', 'Arch', 'Topol', 'Homol']]
# blast_df.columns = ['domain1', 'domain2', 'blast_identity','qcov','E_value', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'Class2', 'Arch2', 'Topol2', 'Homol2']


# blast_df['key_id']=['-'.join(sorted(combine)) for combine in zip(blast_df['domain1'], blast_df['domain2'])]
# tmp=blast_df.groupby('key_id')['blast_identity_ss']
# blast_df['blast_identity_min_ss'] = tmp.transform('min')
# blast_df['blast_identity_max_ss'] = tmp.transform('max')

# tmp=blast_df.groupby('key_id')['qcov']
# blast_df['qcov_min_ss'] = tmp.transform('min')
# blast_df['qcov_max_ss'] = tmp.transform('max')

# tmp=blast_df.groupby('key_id')['E_value']
# blast_df['E_value_min_ss'] = tmp.transform('min')
# blast_df['E_value_max_ss'] = tmp.transform('max')



# blast_df.drop_duplicates(subset='key_id', keep="last", inplace=True)



# blast_df['cath_superFamily'] = np.where(blast_df["Homol1"] == blast_df["Homol2"], 1, 0)



# ######## plots

# #### blast identity min 

# same_homo = blast_df.query('cath_superFamily == 1') ['blast_identity_min']
# diff_homo = blast_df.query('cath_superFamily == 0') ['blast_identity_min']
# same_homo_df = pd.DataFrame({'blast_identity_min':same_homo.values, 'Catagory':'Same Superfamily'})
# diff_homo_df = pd.DataFrame({'blast_identity_min':diff_homo.values, 'Catagory':'Different Superfamily'})
# same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

# g2 = sns.violinplot( y=same_homo_not_same_homo_df["blast_identity_min"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# #plt.ylim(0, 1)
# #g2.set(ylabel=None)
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/blast_identity_min_SS_super_safe_cath_after_ck1-2-3.png', format="png")
# plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/blast_identity_min_SS_super_safe_cath_after_ck1-2-3.svg', format="svg")
# plt.close()
