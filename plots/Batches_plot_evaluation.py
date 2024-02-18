import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, auc,roc_auc_score
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report


print("batch 24")

i=24
df_SS = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/distance_comparison_SS_batch-U_'+str(i)+'.csv')
df_SS['seq1_len']=df_SS['sequence1'].str.len()
df_SS['seq2_len']=df_SS['sequence2'].str.len()

df_SS['key_id']=['-'.join(sorted(combine)) for combine in zip(df_SS['id'], df_SS['dimain2'])]
df_SS=df_SS[['key_id','distance','seq1_len', 'seq2_len']]

df_TM = pd.read_csv('/group/bioinf_protstr/Ballal/TM_merge_all/batch'+str(i)+'/tm_merge_batch'+str(i)+'.csv')
df_TM['key_id']=['-'.join(sorted(combine)) for combine in zip(df_TM['query'], df_TM['subject'])]
df_TM=df_TM[['key_id','tm_query','tm_subject']]

df_SS=df_SS.merge(df_TM, on='key_id', how='inner')

df_SS["tm_min"] = df_SS[["tm_query", "tm_subject"]].min(axis=1)
df_SS['max_seq_len'] = df_SS[["seq1_len", "seq2_len"]].max(axis=1)
df_SS['lev%_SS'] = df_SS['distance']/ df_SS['max_seq_len']
#### x
df_SS['seq_idnt'] = abs(1-df_SS["lev%_SS"])

print("",df_SS.columns.tolist())



# ####### violin plot all for tm min

print("violin plot all")


df_SS["tm_min"] = df_SS["tm_min"].round(decimals = 1)

seq_dist_SS = df_SS['seq_idnt'].to_list()
tm_score_min = df_SS["tm_min"].to_list()



temp_dict = {'tm_score_min':tm_score_min, 'seq_ident_SS':seq_dist_SS}
temp_df = pd.DataFrame(temp_dict)
print("temp_df........len",len(temp_df))

plt.figure(figsize=(20,10))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
sns.violinplot(data=temp_df, x='tm_score_min', y='seq_ident_SS')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_batches_plots/tm_min_lev%_violin_all_original_SS_batch_'+str(i)+'.svg', format="svg")
plt.close()



###  Evaluation


df_SS["lev_class"] = np.where(df_SS["lev%_SS"] >=0.5, 0, 1)
df_SS["tm_min_class"] = np.where(df_SS["tm_min"] >=0.5, 1, 0)
df_SS["seq_idnt_class"] = np.where(df_SS['seq_idnt'] >=0.5, 1, 0)

fpr, tpr, thresholds = roc_curve(df_SS["tm_min_class"],df_SS['seq_idnt'], pos_label=1)
tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2
best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]


plt.plot(thresholds,tpr, color="blue",label="TPR" )
plt.plot(thresholds,fpr, color="red",label="FPR" )
plt.plot(thresholds,tnr, color="green", label="TNR")
plt.plot(thresholds,fnr, color="orange",label="FNR" )
plt.plot(thresholds,balanced_accuracy, color="black",label="BA")
plt.legend(loc="lower right")

plt.axvline(x=best_thre[0], color="darkgrey",linestyle='--')
plt.figtext(best_thre, 0, best_thre, wrap=True, horizontalalignment='center', fontsize=12)

plt.axhline(y=round(max(balanced_accuracy),2), color="darkgrey",linestyle='--')
plt.figtext(0.08, 0.78, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)

plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_batches_plots/tm_min_lev_balanced_accuracy_batch_'+str(i)+'.svg', format="svg")
plt.close()


####for table

# tn, fp, fn, tp = confusion_matrix(df_SS["tm_min_class"],df_SS["seq_idnt_class"],labels=[0,1]).ravel()
# print("tn, fp, fn, tp for seq_idnt 0.5 = ",tn, fp, fn, tp)

# tpr = tp/(tp+fn)
# tnr= tn/(tn+fp)
# fpr = fp/(fp+tn)
# fnr = fn/(fn+tp)
# balanced_accuracy=(tpr+tnr)/2

# print("batch 0",'fpr',fpr,'tpr',tpr,'tnr',tnr,'fnr',fnr ,'balanced accuracy',balanced_accuracy)

# print("batch 0",classification_report(df_SS["tm_min_class"],df_SS["seq_idnt_class"],labels=[0,1]))





