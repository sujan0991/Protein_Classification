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


## cath s20

df_full = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s20_results.csv',usecols=['domain_p2','domain_p1','ss_score','TM_min','cath_superFamily','key_id'])
df_reduced = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_s20_LEV_reduced_results.csv')
df_reduced['key_id']=['-'.join(sorted(combine)) for combine in zip(df_reduced['id'], df_reduced['query_P'])]
df_reduced = df_reduced[['SS_distance', 'key_id']]  

df_full = df_full.merge(df_reduced,left_on='key_id',right_on='key_id')


df_seq=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/cath_non_redundant_pdbs_s20_sequence_reduced_ForkPoolWorker-1.csv')
df_seq = df_seq[~df_seq.AA_seq.isna()]
df2 = df_seq.query('AA_seq.str.contains("\?")')
df2_ids = df2['id'].tolist()
df_seq = df_seq[~df_seq['id'].isin(df2_ids)]
df_seq = df_seq[['id','SS_seq', 'SS_seq_reduced']]
df_seq['ss_reduced_length'] = df_seq['SS_seq_reduced'].str.len().astype(int)


df_full.columns = ['domain_p2', 'id', 'ss_score','key_id','TM_min',  'cath_superFamily',
       'SS_distance']
df_full = pd.merge(df_full, df_seq ,how='left', on="id")
df_full.columns = ['id', 'domain_p1', 'ss_score','key_id', 'TM_min', 'cath_superFamily',
       'SS_distance', 'SS_seq_p1', 'SS_seq_reduced_p1', 'ss_reduced_length_p1']
df_full = pd.merge(df_full, df_seq ,how='left', on="id")

df_full.columns = ['domain_p2', 'domain_p1', 'ss_score','key_id', 'TM_min',  'cath_superFamily',
       'SS_distance', 'SS_seq_p1', 'SS_seq_reduced_p1', 'ss_reduced_length_p1','SS_seq_p2', 
       'SS_seq_reduced_p2', 'ss_reduced_length_p2']



df_full['max_seq_len'] = df_full[["ss_reduced_length_p1", "ss_reduced_length_p2"]].max(axis=1)
df_full['levIdent_SS'] = pd.to_numeric(df_full['SS_distance'])/ df_full['max_seq_len']
df_full['ss_score_reduced'] = abs(1-df_full["levIdent_SS"]).astype(np.float16)

df_full = df_full[['domain_p2', 'domain_p1', 'ss_score', 'ss_score_reduced','key_id', 'TM_min', 
       'cath_superFamily', 'SS_seq_p1', 'SS_seq_reduced_p1','SS_seq_p2', 'SS_seq_reduced_p2']]
df_full.columns = ['domain_p2', 'domain_p1', 'ss_score_full', 'ss_score_reduced','key_id', 'TM_min', 
       'cath_superFamily', 'SS_seq_full_p1', 'SS_seq_reduced_p1','SS_seq_full_p2', 'SS_seq_reduced_p2']






## query_id,id,local_alignments_full,local_ident_full,query_aln_string_full,target_aln_string_full,local_alignments_reduced,local_ident_reduced,query_aln_string_reduced,target_aln_string_reduced

df_local = pd.read_csv('/group/bioinf/Ballal/cath_s20_local_results_full_reduced.csv',usecols=['query_id','id','local_ident_full','local_ident_reduced'])
df_local['key_id']=['-'.join(sorted(combine)) for combine in zip(df_local['id'], df_local['query_id'])]
df_local = df_local[['local_ident_full','local_ident_reduced', 'key_id']]  


df_full = df_full.merge(df_local,left_on='key_id',right_on='key_id')

df_full = df_full [['domain_p1','domain_p2', 'ss_score_full', 'ss_score_reduced','local_ident_full','local_ident_reduced',
                    'SS_seq_full_p1','SS_seq_full_p2','SS_seq_reduced_p1','SS_seq_reduced_p2','key_id',
                    'TM_min', 'cath_superFamily']]


df_full.to_csv('/group/bioinf/Ballal/cath_s20_ss_local_results_full_reduced.csv', index= False)




## AUC

## full / reduced

fpr, tpr, thresholds = roc_curve(df_full["cath_superFamily"],df_full['TM_min'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...TM_min.",roc_auc)  ### 0.89816031807

plt.plot(fpr,tpr,color="blue",label="TM")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')


# fpr, tpr, thresholds = roc_curve(df_full["cath_superFamily"],df_full['ss_score_reduced'], pos_label=1)
# roc_auc = auc(fpr, tpr)
# print("roc_auc...ss.",roc_auc,type(roc_auc))##full: 0.807628386 reduced : 0.78
# plt.plot(fpr,tpr,color="green",label="SS")

fpr, tpr, thresholds = roc_curve(df_full["cath_superFamily"],df_full['local_ident_reduced'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc....local_ident_reduced",roc_auc)## roc_auc....local_ident_full 0.72 , reduced: 0.670

plt.plot(fpr,tpr,color="orange",label="Local")



plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/tm_min_lev_roc_cath_s20.png', format="png")
plt.savefig('/group/bioinf/Ballal/plots/tm_min_local_reduced_roc_cath_s20.svg', format="svg")
plt.close()










## scop

## protein1,protein2,is_same_scope_Superfamily,SS_score_DSSP,key_id,tm_min,tm_max,AA_score,SS_score,scope_Family
df_full = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_DSSP_LEV_TM_results.csv',usecols=['protein1','protein2','is_same_scope_Superfamily','key_id','SS_score','tm_min'])
df_reduced = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOP_LEV_reduced_results.csv')
df_reduced['key_id']=['-'.join(sorted(combine)) for combine in zip(df_reduced['id'], df_reduced['query_P'])]
df_reduced = df_reduced[['SS_distance', 'key_id']]  
df_full = df_full.merge(df_reduced,left_on='key_id',right_on='key_id')
## 

df_seq=pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_seq_reduced_ForkPoolWorker-1.csv')
df_seq['ss_reduced_length'] = df_seq['SS_seq_reduced'].str.len().astype(int)
df_seq =df_seq [['id','SS_seq', 'SS_seq_reduced','ss_reduced_length']]

##
df_full.columns = ['id', 'protein2', 'is_same_scope_Superfamily', 'key_id', 'tm_min',
       'ss_score', 'SS_distance']
df_full = pd.merge(df_full, df_seq ,how='left', on="id")
df_full.columns = ['protein1', 'id', 'is_same_scope_Superfamily', 'key_id','tm_min',
       'ss_score_full', 'SS_distance', 'SS_seq_p1', 'SS_seq_reduced_p1',
       'ss_reduced_length_p1']
df_full = pd.merge(df_full, df_seq ,how='left', on="id")
df_full.columns = ['protein1', 'protein2', 'is_same_scope_Superfamily', 'key_id', 'tm_min',
       'ss_score_full', 'SS_distance', 'SS_seq_p1', 'SS_seq_reduced_p1',
       'ss_reduced_length_p1', 'SS_seq_p2', 'SS_seq_reduced_p2',
       'ss_reduced_length_2']
df_full['max_seq_len'] = df_full[["ss_reduced_length_p1", "ss_reduced_length_2"]].max(axis=1)
df_full['levIdent_SS'] = pd.to_numeric(df_full['SS_distance'])/ df_full['max_seq_len']
df_full['ss_score_reduced'] = abs(1-df_full["levIdent_SS"]).astype(np.float16)

df_full = df_full[['protein1', 'protein2', 'ss_score_full','ss_score_reduced', 'SS_seq_p1','SS_seq_p2',
                   'SS_seq_reduced_p1','SS_seq_reduced_p2','key_id','tm_min','is_same_scope_Superfamily']]


## scop_local_results_full_reduced.csv   ## have less result, find out why

df_local = pd.read_csv('/group/bioinf/Ballal/scop_local_results_full_reduced.csv',usecols=['query_id','id','local_ident_full','local_ident_reduced'])
df_local['key_id']=['-'.join(sorted(combine)) for combine in zip(df_local['id'], df_local['query_id'])]
df_local = df_local[['local_ident_full','local_ident_reduced', 'key_id']]  




df_local = df_local.merge(df_full,left_on='key_id',right_on='key_id')

df_local = df_local [['protein1','protein2', 'ss_score_full', 'ss_score_reduced','local_ident_full', 'local_ident_reduced',
                       'SS_seq_p1','SS_seq_p2', 'SS_seq_reduced_p1', 'SS_seq_reduced_p2','key_id',
                       'tm_min','is_same_scope_Superfamily']]
df_local.columns = ['protein1','protein2', 'ss_score_full', 'ss_score_reduced','local_ident_full', 'local_ident_reduced',
                       'SS_seq_full_p1','SS_seq_full_p2', 'SS_seq_reduced_p1', 'SS_seq_reduced_p2','key_id',
                       'tm_min','is_same_scope_Superfamily']



df_local.to_csv('/group/bioinf/Ballal/scop_ss_local_results_full_reduced.csv', index= False)


## AUC

## full / reduced

fpr, tpr, thresholds = roc_curve(df_local["is_same_scope_Superfamily"],df_local['tm_min'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc...TM_min.",roc_auc)  ### 

plt.plot(fpr,tpr,color="blue",label="TM")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')


# fpr, tpr, thresholds = roc_curve(df_local["is_same_scope_Superfamily"],df_local['ss_score_reduced'], pos_label=1)
# roc_auc = auc(fpr, tpr)
# print("roc_auc...ss.",roc_auc,type(roc_auc))##full:0.888  reduced : 0.840
# plt.plot(fpr,tpr,color="green",label="SS")

fpr, tpr, thresholds = roc_curve(df_local["is_same_scope_Superfamily"],df_local['local_ident_reduced'], pos_label=1)
roc_auc = auc(fpr, tpr)
print("roc_auc....local_ident_reduced",roc_auc)## roc_auc....local_ident_full 0.768  , reduced: 0.661

plt.plot(fpr,tpr,color="orange",label="Local")



plt.axline((1, 1), slope=1,color="darkgrey",linestyle='--')

plt.legend(loc="lower right")
#plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/tm_min_lev_roc_cath_s20.png', format="png")
plt.savefig('/group/bioinf/Ballal/plots/tm_min_local_reduced_roc_scop.svg', format="svg")
plt.close()












### manual testing in scop


from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Levenshtein import distance as lev


def perform_alignment(query_string, target_string, c=[2, -1, -3, -2]):
    shorter_seq_length = min(len(query_string),len(target_string))
    local_alignments = pairwise2.align.localms(query_string, target_string, c[0], c[1], c[2], c[3], one_alignment_only=True)
    if len(local_alignments) > 0:local_alignments = local_alignments[0][2]
    else: local_alignments =0
    local_ident = local_alignments / (shorter_seq_length*c[0])
    ## lev 
    longer_seq_length = max(len(query_string),len(target_string))
    lev_score = lev(query_string, target_string)
    ss_score = abs(1 - lev_score / longer_seq_length)
    return ss_score, local_ident


query_string = 'LLLLLHHHHHHHHHHHHHHLLLLLLLSSSSSSLLLLLLLLLLSSSSSSSLLLLSSSSSSSLLLHHHHHHHHHHHHHHHLLLLLLL'
target_string = 'LHHHHHHHHHHHHHHHHHHHHHHHLSSSLLLLLSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSLLLLLHHHHHHHHHHLLLLSSSSSSSSSSSLLLLLLLLLLLLLSSSSSSSSSSLLLLLLLHHHHHHHHHHHHHHHHHHHHHHHHHLLLLLLLLLLLSSSSHHHHHHHLLLLLHHHHHHHHHHHHLSSSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSLLLLSSSSSSSSSSLLLHHHHHHHHHHHLLHHHHHLHHHHHHHHLLLLLSSSSSSSHHHHHHHHHLLLLHHHHLLLLLLHHHHHHLLLLL'
perform_alignment(query_string, target_string)


