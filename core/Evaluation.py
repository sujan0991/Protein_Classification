import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc,roc_auc_score
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report


lev_tm_df = pd.read_csv('lev_tm_10k_original_SS_df.csv')

print("lev_tm_df len",len(lev_tm_df))


# excludeList=pd.read_csv('/group/bioinf_protstr/Ballal/excludeList.txt',header=None)
# >>> len(lev_tm_df)
# 11255140

# excludeList=list(excludeList[0])
# >>> lev_tm_df = lev_tm_df[~(lev_tm_df.domain1_x.isin(excludeList)) & ~(lev_tm_df.domain2_x.isin(excludeList))]
# >>> len(lev_tm_df)
# 8555316


lev_tm_df['seq_1_len'] = lev_tm_df['sequence1'].apply(len)
lev_tm_df['seq_2_len'] = lev_tm_df['sequence2'].apply(len)
lev_tm_df['max_seq_len'] = lev_tm_df[["seq_1_len", "seq_2_len"]].max(axis=1)
lev_tm_df['lev%_SS'] = lev_tm_df['lev_distance']/ lev_tm_df['max_seq_len']
lev_tm_df['lev%_AA'] = lev_tm_df['lev_distance_AA']/ lev_tm_df['max_seq_len']


#### x
lev_tm_df['seq_idnt_SS'] = abs(1-lev_tm_df["lev%_SS"])
lev_tm_df['seq_idnt_AA'] = abs(1-lev_tm_df["lev%_AA"])
lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)

lev_tm_df = lev_tm_df[['domains','Homol1', 'Homol2','tm_min','seq_idnt_SS','seq_idnt_AA']]
###### set everything to 0 and 1

##lev_tm_df["lev_class"] = np.where(lev_tm_df["lev%_SS"] >=0.5, 0, 1)
lev_tm_df["tm_min_class"] = np.where(lev_tm_df["tm_min"] >=0.5, 1, 0)
lev_tm_df["seq_idnt_class_SS"] = np.where(lev_tm_df['seq_idnt'] >=0.5, 1, 0)
lev_tm_df["seq_idnt_class_AA"] = np.where(lev_tm_df['seq_idnt_AA'] >=0.3, 1, 0)
lev_tm_df["cath_superFamily"] = np.where(lev_tm_df["Homol1"] == lev_tm_df["Homol2"], 1, 0)


### confusion_matrix

# tn, fp, fn, tp = confusion_matrix(lev_tm_df["tm_min_class"], lev_tm_df["lev_class"],)
# print("tn, fp, fn, tp for lev% 0.5 = ",tn, fp, fn, tp)




fpr, tpr, thresholds = roc_curve(lev_tm_df["tm_min_class"],lev_tm_df['tm_min'], pos_label=1)
print("fpr, tpr len",len(fpr),len(tpr))
roc_auc = auc(fpr, tpr)
print("roc_auc....",roc_auc)


# #create ROC curve
plt.plot(fpr,tpr,color="red",label="TM")
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc="lower right")
plt.savefig('tm_min_lev_roc_test_set.png', format="png")
plt.savefig('tm_min_lev_roc_test_set.svg', format="svg")
plt.close()

##fpr, tpr, thresholds = roc_curve(df_SS["tm_min_class"],df_SS['seq_idnt'], pos_label=1)

#fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['seq_idnt'], pos_label=1)
#fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['seq_idnt_AA'], pos_label=1)
fpr, tpr, thresholds = roc_curve(lev_tm_df["cath_superFamily"],lev_tm_df['tm_min'], pos_label=1)

tnr=1-fpr
fnr=1-tpr
balanced_accuracy=(tpr+tnr)/2

best_thre=thresholds[balanced_accuracy==max(balanced_accuracy)]
balanced_accuracy_05= balanced_accuracy[thresholds==0.5]
balanced_accuracy_03= balanced_accuracy[thresholds==0.3]

plt.axis([0, 1, 0, 1])
plt.plot(thresholds,tpr, color="blue",label="tpr" )
plt.plot(thresholds,fpr, color="red",label="fpr" )
plt.plot(thresholds,tnr, color="green", label="tnr")
plt.plot(thresholds,fnr, color="orange",label="fnr" )
plt.plot(thresholds,balanced_accuracy, color="black",label="ba")
plt.legend(loc="lower right")

plt.axvline(x=best_thre[0], color="darkgrey",linestyle='--')
plt.figtext(0.27, 0.03, round(best_thre[0],2), wrap=True, horizontalalignment='center', fontsize=12)

plt.axhline(y=round(max(balanced_accuracy),2), color="darkgrey",linestyle='--')
plt.figtext(0.08, 0.78, round(max(balanced_accuracy),2), wrap=True, horizontalalignment='center', fontsize=12)


plt.axvline(x=thresholds[(thresholds>=0.5)&(thresholds<=0.5)], color="darkgrey",linestyle='--')

plt.axhline(y=round(balanced_accuracy_05[0],2), color="darkgrey",linestyle='--')
plt.figtext(0.08, 0.74, round(balanced_accuracy_05[0],2), wrap=True, horizontalalignment='center', fontsize=12)

# plt.axhline(y=round(balanced_accuracy_03[0],2), color="darkgrey",linestyle='--')
# plt.figtext(0.08, 0.65, round(balanced_accuracy_03[0],2), wrap=True, horizontalalignment='center', fontsize=12)



# plt.savefig('tm_min_superfamily_balanced_accuracy_test_set_2.png', format="png")
plt.savefig('tm_min_superfamily_balanced_accuracy_test_set_2.png', format="png")
plt.savefig('tm_min_superfamily_balanced_accuracy_test_set_2.svg', format="svg")

plt.close()


####for table

tn, fp, fn, tp = confusion_matrix(lev_tm_df["tm_min_class"],lev_tm_df["seq_idnt_class_SS"],labels=[0,1]).ravel()
print("tn, fp, fn, tp for seq_idnt 0.5 = ",tn, fp, fn, tp)

tpr = tp/(tp+fn)
tnr= tn/(tn+fp)
fpr = fp/(fp+tn)
fnr = fn/(fn+tp)
balanced_accuracy=(tpr+tnr)/2

print("batch 0",'fpr',fpr,'tpr',tpr,'tnr',tnr,'fnr',fnr ,'balanced accuracy',balanced_accuracy)

print("batch 0",classification_report(lev_tm_df["tm_min_class"],lev_tm_df["seq_idnt_class_SS"],labels=[0,1]))


exit()




lev_tm_df["seq_idnt_class"] = np.where(x >=0.5, 1, 0)

print(lev_tm_df.seq_idnt_class.unique())

tn, fp, fn, tp = confusion_matrix(lev_tm_df["tm_min_class"],lev_tm_df["seq_idnt_class"],labels=[0,1]).ravel()
print("tn, fp, fn, tp for seq_idnt 0.5 = ",tn, fp, fn, tp)

tpr = tp/(tp+fn)
tnr= tn/(tp+fn)
fpr = fp/(fp+tn)
fnr = fn/(fn+tp) 


balanced_accuracy=(tpr+tnr)/2

print('fpr',fpr,'tpr',tpr,'tnr',tnr,'fnr',fnr ,'balanced accuracy',balanced_accuracy)

print(classification_report(lev_tm_df["tm_min_class"],lev_tm_df["seq_idnt_class"],labels=[0,1]))
