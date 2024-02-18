import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

#define function to swap columns
def swap_columns(df, col1, col2):
    col_list = list(df.columns)
    x, y = col_list.index(col1), col_list.index(col2)
    col_list[y], col_list[x] = col_list[x], col_list[y]
    df = df[col_list]
    return df



########### lev tm merge

lev_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/distance_comparison_test_set_2_10k_superfamily_domain_original_SS.csv')
lev_original_AA_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/distance_comparison_test_set_2_10k_without_AA_original_AA.csv') 
print("",lev_original_AA_df.columns.tolist())
tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/tm_merge_10k_test_set_2.csv')


lev_df.columns = ['index', 'domain1','sequence1', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'lev_distance', 'domain2','sequence2', 'Class2', 'Arch2', 'Topol2', 'Homol2']
lev_original_AA_df.columns = ['index', 'domain1', 'domain2', 'lev_distance_AA']
tm_df.columns = ['domain1', 'domain2', 'length', 'rmsd', 'tm_query', 'tm_subject', 'cigar']


lev_df['domains'] = ['-'.join(sorted(combine)) for combine in zip(lev_df['domain1'], lev_df['domain2'])]
lev_original_AA_df['domains'] = ['-'.join(sorted(combine)) for combine in zip(lev_original_AA_df['domain1'], lev_original_AA_df['domain2'])]
tm_df['domains'] = ['-'.join(sorted(combine)) for combine in zip(tm_df['domain1'], tm_df['domain2'])]


lev_tm_df = lev_df.merge(tm_df,left_on='domains',right_on='domains')
lev_tm_df = lev_tm_df.merge(lev_original_AA_df,left_on='domains',right_on='domains')



# #print("///////",lev_tm_df.columns.tolist())

# lev_tm_df = lev_tm_df[lev_tm_df['lev_distance'].notna()]
# lev_tm_df = lev_tm_df[lev_tm_df['tm_query'].notna()]
# lev_tm_df = lev_tm_df[lev_tm_df['lev_distance_AA'].notna()]

# print("lev_tm_df len",len(lev_tm_df))
lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)

lev_tm_df['seq1_len']=lev_tm_df['sequence1'].str.len()
lev_tm_df['seq2_len']=lev_tm_df['sequence2'].str.len()
lev_tm_df['max_seq_len'] = lev_tm_df[["seq1_len", "seq2_len"]].max(axis=1)
# df_SS['min_seq_len'] = df_SS[["seq1_len", "seq2_len"]].min(axis=1)
 
lev_tm_df['lev%_SS'] = lev_tm_df['lev_distance']/ lev_tm_df['max_seq_len']
lev_tm_df['lev%_AA'] = lev_tm_df['lev_distance_AA']/ lev_tm_df['max_seq_len']

#### x
lev_tm_df['seq_idnt'] = abs(1-lev_tm_df["lev%_SS"])
lev_tm_df['seq_idnt_AA'] = abs(1-lev_tm_df["lev%_AA"])



lev_tm_df = lev_tm_df[['domains','Homol1', 'Homol2','tm_min','seq_idnt','seq_idnt_AA']]

lev_tm_df["cath_superFamily"] = np.where(lev_tm_df["Homol1"] == lev_tm_df["Homol2"], 1, 0)
lev_tm_df = lev_tm_df[['domains','cath_superFamily','tm_min','seq_idnt','seq_idnt_AA']]

# lev_tm_df.to_csv('lev_tm_10k_original_SS_df.csv')



exit()



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




## df2 = pd.read_csv("cath-domain-description-v4_3_0.csv")


df_SS=df_SS.merge(df_TM, on='key_id', how='inner')


# lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_tm_10k_original_SS_df.csv')

df_SS["tm_min"] = df_SS[["tm_query", "tm_subject"]].min(axis=1)
# df_SS["tm_max"] = df_SS[["tm_query", "tm_subject"]].max(axis=1)

print("",df_SS.columns.tolist())

print("lev_tm_df len",len(df_SS))



#### seq identity plot ## see in density plot section

# lev_tm_df['seq1_len'] = lev_tm_df['sequence1'].apply(len)
 
# lev_tm_df['seq2_len'] = lev_tm_df['sequence2'].apply(len)

# lev_tm_df['seq_len_diff'] = abs(lev_tm_df['seq1_len'] - lev_tm_df['seq2_len']) 

df_SS['max_seq_len'] = df_SS[["seq1_len", "seq2_len"]].max(axis=1)
# df_SS['min_seq_len'] = df_SS[["seq1_len", "seq2_len"]].min(axis=1)
 
df_SS['lev%_SS'] = df_SS['distance']/ df_SS['max_seq_len']

# df_SS['lev%_AA'] = df_SS['lev_distance_AA']/ df_SS['max_seq_len']

#### x
df_SS['seq_idnt'] = abs(1-df_SS["lev%_SS"])
# df_SS['seq_idnt_AA'] = abs(1-df_SS["seq_idnt_AA"])

print("",df_SS.columns.tolist())



# lev_dist_SS = df_SS['lev_distance'].to_list()
# lev_dist_AA = df_SS['lev_distance_AA'].to_list()
# seq_dist_SS = df_SS['lev%_SS'].to_list()
# seq_dist_AA = df_SS['lev%_AA'].to_list()
# tm_score_min = df_SS["tm_min"].to_list()
# tm_score_max = df_SS["tm_max"].to_list()

##print(lev_tm_df.isnull().any()) ########## check if any column have nan value

# lev_dist_SS = np.nan_to_num(lev_dist_SS)
# lev_dist_AA = np.nan_to_num(lev_dist_AA)
# seq_dist_SS = np.nan_to_num(seq_dist_SS)
# tm_score_min = np.nan_to_num(tm_score_min)
# seq_dist_AA = np.nan_to_num(seq_dist_AA)
# tm_score_max = np.nan_to_num(tm_score_max)





# ####### violin plot all for tm min

print("violin plot all")


lev_tm_df["tm_min"] = lev_tm_df["tm_min"].round(decimals = 1)
# lev_tm_df["tm_max"] = lev_tm_df["tm_max"].round(decimals = 1)


seq_dist_SS = lev_tm_df['seq_idnt'].to_list()
# seq_dist_AA = df_SS['seq_idnt_AA'].to_list()
tm_score_min = lev_tm_df["tm_min"].to_list()
# tm_score_max = df_SS["tm_max"].to_list()


tm_0 = 0
tm_1 =0
tm_2 = 0
tm_3 = 0
tm_4 =0
tm_5 = 0
tm_6 = 0
tm_7 =0
tm_8 = 0
tm_9 = 0
tm_10 =0


for i,value in enumerate(tm_score_min):
    if tm_score_min[i] == 0.0:
        tm_0 = tm_0 + 1
    elif tm_score_min[i] == 0.1:
        tm_1 = tm_1 + 1 
    elif tm_score_min[i] == 0.2:
        tm_2 = tm_2 + 1
    elif tm_score_min[i] == 0.3:
        tm_3 = tm_3 + 1
    elif tm_score_min[i] == 0.4:
        tm_4 = tm_4 + 1
    elif tm_score_min[i] == 0.5:
        tm_5 = tm_5 + 1
    elif tm_score_min[i] == 0.6:
        tm_6 = tm_6 + 1
    elif tm_score_min[i] == 0.7:
        tm_7 = tm_7 + 1
    elif tm_score_min[i] == 0.8:
        tm_8 = tm_8 + 1  
    elif tm_score_min[i] == 0.9:
        tm_9 = tm_9 + 1 
    else :                               
        tm_10 = tm_10 + 1
        
print("tm-min 0-1",tm_0,tm_1,tm_2,tm_3,tm_4,tm_5,tm_6,tm_7,tm_8,tm_9,tm_10) 




temp_dict = {'tm_score_min':tm_score_min, 'seq_ident_SS':seq_dist_SS}

temp_df = pd.DataFrame(temp_dict)

print("temp_df........len",len(temp_df))

plt.figure(figsize=(20,8))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
sns.violinplot(data=temp_df, x='tm_score_min', y='seq_ident_SS')
plt.savefig('tm_min_lev%_violin_all_original_SS.png', format="png")
plt.savefig('tm_min_lev%_violin_all_original_SS.svg', format="svg")
plt.close()


########### for AA sequence

# temp_dict = {'tm_score_min':tm_score_min, 'seq_ident_AA':seq_dist_AA}

# temp_df = pd.DataFrame(temp_dict)

# print("temp_df........len",len(temp_df))
# plt.figure(figsize=(20,8))
# sns.violinplot(data=temp_df, x='tm_score_min', y='seq_ident_AA')
# plt.savefig('tm_min_lev%_violin_all_original_AA_test_set_2.png', format="png")
# plt.close()






############### correlations ##################

# corr, p_v = pearsonr(lev_dist_SS,lev_dist_AA)
# print('Pearsons correlation lev_dist_AA vs lev_dist_SS all: %.3f' % corr, p_v)

# corr, p_v = pearsonr(seq_dist_SS,seq_dist_AA)
# print('Pearsons correlation seq_dist_SS vs seq_dist_AA all: %.3f' % corr, p_v)


# corr, p_v = pearsonr(seq_dist_SS,tm_score_min)
# print('Pearsons correlation TM-min vs lev_SS all: %.3f' % corr, p_v)

# corr, p_v = pearsonr(seq_dist_AA,tm_score_min)
# print('Pearsons correlation TM-min vs lev_AA all: %.3f' % corr, p_v)


# corr, p_v = pearsonr(seq_dist_SS,tm_score_max)
# print('Pearsons correlation TM-max vs lev_SS all: %.3f' % corr, p_v)


# corr, p_v = pearsonr(seq_dist_AA,tm_score_max)
# print('Pearsons correlation TM-max vs lev_AA all: %.3f' % corr, p_v)

corr, p_v = pearsonr(seq_dist_SS[tm_score_min>=0.1],tm_score_min[tm_score_min>=0.1])
print('Pearsons correlation TM-min vs lev_SS all: %.3f' % corr, p_v)





# # # ########### lev superfamily for SS/AA

##'domain1','domain2','AA_distance', 'SS_distance', 'domain1_length', 'Class1', 'Arch1', 'Topol1', 'Homol1','domain2_length', 'Class2', 'Arch2', 'Topol2', 'Homol2'

df_merged['max_seq_len'] = df_merged[["domain1_length", "domain2_length"]].max(axis=1)
df_merged['lev%_SS'] = df_merged['SS_distance']/ df_merged['max_seq_len']
#### x
df_merged['seq_idnt'] = abs(1-df_merged["lev%_SS"])

same_homo = df_merged.query('Homol1 == Homol2') ['seq_idnt']
diff_homo = df_merged.query('Homol1 != Homol2') ['seq_idnt']

below_point_5 = 0
above_point_5 = 0
for i, v in same_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1   

print('lev for same superfamily:total',len(same_homo), 'above point 0.5', above_point_5 ,'below point 0.5', below_point_5)


below_point_5 = 0
above_point_5 = 0
for i, v in diff_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('lev for different superfamily:total',len(diff_homo), ' above point 0.5', above_point_5 ,'below point 0.5', below_point_5)
 
same_homo_df = pd.DataFrame({'seq_idnt_SS':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'seq_idnt_SS':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["seq_idnt_SS"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#g2.set(ylabel=None)
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/seq_idnt_SS_super_safe_cath.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/seq_idnt_SS_super_safe_cath.svg', format="svg")
plt.close()



# ############ Same and different Superfamily for tm max/tm min

##df_SS.loc[(df_SS['Homol1']==df_SS['Homol2']) & (df_SS['tm_min']< 0.4),['domains','tm_min','Homol1','Homol2']]

same_homo = lev_tm_df.query('Homol1 == Homol2') ['tm_min']
diff_homo = lev_tm_df.query('Homol1 != Homol2') ['tm_min']

below_point_5 = 0
above_point_5 = 0

for i, v in same_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1    
print('tm for same superfamily:total',len(same_homo), 'above point 0.5', above_point_5 ,'below point =0.5', below_point_5)


below_point_5 = 0
above_point_5 = 0

for i, v in diff_homo.items():
    if v > 0.5:
        above_point_5 = above_point_5 + 1
    else:
        below_point_5 = below_point_5 + 1         

print('tm for different superfamily:total',len(diff_homo), ' above point 0.5', above_point_5 ,'below point =0.5', below_point_5)
 

same_homo_df = pd.DataFrame({'TM_min':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'TM_min':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["TM_min"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#g2.set(ylabel=None)
plt.savefig('tm_min_super_10k_original_AA_test_set_2_v.png', format="png")
plt.close()





##############  tm < 0.5 and tm > 0.5 for SS/AA

group1 = lev_tm_df.query('tm_min < 0.5') ['seq_idnt']
group2 = lev_tm_df.query('tm_min >= 0.5') ['seq_idnt']

print('tm vs seq_ident_SS:total',len(group1)+len(group2), ' above point 0.5', len(group2) ,'below point =0.5',len(group1))
 
group1_df = pd.DataFrame({'seq_ident_SS':group2.values, 'Catagory':'tm_min >= 0.5'})
group2_df = pd.DataFrame({'seq_ident_SS':group1.values, 'Catagory':'tm_min < 0.5'})

group1_group2_df = group1_df.append(group2_df, ignore_index=True)

g2 = sns.violinplot( y=group1_group2_df["seq_ident_SS"], x=group1_group2_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.savefig('tm_lev_violin_below_above_5_original_SS_test_set_2.png', format="png")
plt.savefig('tm_lev_violin_below_above_5_original_SS_test_set_2.svg', format="svg")
plt.close() 






# #### plot againest max length

# max_seq_len = lev_tm_df['max_seq_len'].to_list()

# temp_dict = {'tm_score_min':tm_score_min, 'max_seq_len':max_seq_len}

# temp_df = pd.DataFrame(temp_dict)

# print("temp_df........len",len(temp_df))

# plt.figure(figsize=(20,8))
# sns.violinplot(data=temp_df, x='tm_score_min', y='max_seq_len')
# plt.savefig('tm_min_max_seq_len_violin_all_original_SS_test_set_2.png', format="png")
# plt.close()


# ############# againest min seq length

# ####

# min_seq_len = lev_tm_df['min_seq_len'].to_list()

# temp_dict = {'tm_score_min':tm_score_min, 'min_seq_len':min_seq_len}

# temp_df = pd.DataFrame(temp_dict)

# print("temp_df........len",len(temp_df))

# plt.figure(figsize=(20,8))
# sns.violinplot(data=temp_df, x='tm_score_min', y='min_seq_len')
# plt.savefig('tm_min_min_seq_len_violin_all_original_SS_test_set_2.png', format="png")
# plt.close()





# #### plot againest difference in length

# seq_len_diff = lev_tm_df['seq_len_diff'].to_list()

# temp_dict = {'tm_score_min':tm_score_min, 'seq_len_diff':seq_len_diff}

# temp_df = pd.DataFrame(temp_dict)

# print("temp_df........len",len(temp_df))

# plt.figure(figsize=(20,8))
# sns.violinplot(data=temp_df, x='tm_score_min', y='seq_len_diff')
# plt.savefig('tm_min_seq_len_diff_violin_all_original_SS_test_set_2.png', format="png")
# plt.close()





############# againest min seq length

min_seq_len = lev_tm_df['min_seq_len'].to_list()

lev_SS = lev_tm_df['lev%_SS'].to_list()

temp_dict = {'tm_score_min':tm_score_min, 'min_seq_len':min_seq_len, 'lev_SS':lev_SS}

temp_df = pd.DataFrame(temp_dict)

print("temp_df........len",len(temp_df))

plt.figure(figsize=(20,8))
sns.violinplot(data=temp_df[temp_df.min_seq_len>100], x='tm_score_min', y='lev_SS')
plt.savefig('tm_min_min_seq_len_lev_violin_all_original_SS_test_set_2.png', format="png")
plt.close()





####### when tm score 0.0


tm_0 = lev_tm_df[lev_tm_df.tm_min == 0.0]


same_homo = tm_0.query('Homol1 == Homol2') ['lev%_SS']
diff_homo = tm_0.query('Homol1 != Homol2') ['lev%_SS']


print("same superfaamily,diff superfamily",len(same_homo),len(diff_homo))

same_homo_df = pd.DataFrame({'seq_ident_SS':same_homo.values, 'Catagory':'Same Superfamily'})
diff_homo_df = pd.DataFrame({'seq_ident_SS':diff_homo.values, 'Catagory':'Different Superfamily'})
 
same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

g2 = sns.violinplot( y=same_homo_not_same_homo_df["seq_ident_SS"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#g2.set(ylabel=None)
plt.savefig('tm_0_vs_seq_idnt_10k_original_SS_test_set_2_v.png', format="png")
plt.close()





exit()

############ for tm max

tm_0 = 0
tm_1 =0
tm_2 = 0
tm_3 = 0
tm_4 =0
tm_5 = 0
tm_6 = 0
tm_7 =0
tm_8 = 0
tm_9 = 0
tm_10 =0


for i,value in enumerate(tm_score_min):

    if tm_score_min[i] == 0.0:
        tm_0 = tm_0 + 1
    elif tm_score_min[i] == 0.1:
        tm_1 = tm_1 + 1 
    elif tm_score_min[i] == 0.2:
        tm_2 = tm_2 + 1
        
    elif tm_score_min[i] == 0.3:
        tm_3 = tm_3 + 1
        
    elif tm_score_min[i] == 0.4:
        tm_4 = tm_4 + 1
    elif tm_score_min[i] == 0.5:
        tm_5 = tm_5 + 1
    elif tm_score_min[i] == 0.6:
        tm_6 = tm_6 + 1
    elif tm_score_min[i] == 0.7:
        tm_7 = tm_7 + 1
    elif tm_score_min[i] == 0.8:
        tm_8 = tm_8 + 1  
    elif tm_score_min[i] == 0.9:
        tm_9 = tm_9 + 1 
    else :                               
        tm_10 = tm_10 + 1
        
print("tm-max 0-1",tm_0,tm_1,tm_2,tm_3,tm_4,tm_5,tm_6,tm_7,tm_8,tm_9,tm_10) 


temp_dict = {'tm_score_max':tm_score_max, 'seq_ident_SS':seq_dist_SS}

temp_df = pd.DataFrame(temp_dict)

print("temp_df........len",len(temp_df))

sns.violinplot(data=temp_df, x='tm_score_max', y='seq_ident_SS')
plt.savefig('tm_max_lev%_violin_all_original_SS_test_set_2.png', format="png")
plt.close()


temp_dict = {'tm_score_max':tm_score_max, 'seq_ident_AA':seq_dist_AA}

temp_df = pd.DataFrame(temp_dict)

print("temp_df........len",len(temp_df))

sns.violinplot(data=temp_df, x='tm_score_max', y='seq_ident_AA')
plt.savefig('tm_max_lev%_violin_all_original_AA_test_set_2.png', format="png")
plt.close()


exit()



####### when tm score 0.0

# lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)

# lev_tm_df["tm_min"] = lev_tm_df["tm_min"].round(decimals = 1)

# tm_0 = lev_tm_df[lev_tm_df.tm_min == 0.0]


# same_homo = tm_0.query('Homol1 == Homol2') ['lev_distance']
# diff_homo = tm_0.query('Homol1 != Homol2') ['lev_distance']


# print("same superfaamily,diff superfamily",len(same_homo),len(diff_homo))

# same_homo_df = pd.DataFrame({'lev_distance':same_homo.values, 'Catagory':'Same Superfamily'})
# diff_homo_df = pd.DataFrame({'lev_distance':diff_homo.values, 'Catagory':'Different Superfamily'})
 
# same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

# g2 = sns.violinplot( y=same_homo_not_same_homo_df["lev_distance"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# #g2.set(ylabel=None)
# plt.savefig('tm_0_10k_without_AA_test_set_2_v.png', format="png")
# plt.close()






########## 

print("???????????")
group1 = lev_tm_df.query('tm_min >= 0.5 & Homol1 == Homol2') ['lev_distance']
group2 = lev_tm_df.query('tm_min < 0.5 & Homol1 == Homol2') ['lev_distance']
group3 = lev_tm_df.query('tm_min >= 0.5 & Homol1 != Homol2') ['lev_distance']
group4 = lev_tm_df.query('tm_min < 0.5 & Homol1 != Homol2') ['lev_distance']


group1_df = pd.DataFrame({'lev_distance':group1.values, 'Catagory':'tm >= 0.5 & same S.F'})
group2_df = pd.DataFrame({'lev_distance':group2.values, 'Catagory':'tm < 0.5 & same S.F'})
group3_df = pd.DataFrame({'lev_distance':group3.values, 'Catagory':'tm >= 0.5 & diff. S.F'})
group4_df = pd.DataFrame({'lev_distance':group4.values, 'Catagory':'tm < 0.5 & diff. S.F'})

print(len(group1_df),len(group2_df),len(group3_df),len(group4_df))

group1_group2_df = pd.concat([group1_df,group2_df,group3_df,group4_df])


g2 = sns.violinplot( y=group1_group2_df["lev_distance"], x=group1_group2_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
plt.savefig('tm_lev_violin_below_above_5_same_diff_super_original_SS_test_set_2.png', format="png")
plt.close() 

print("???????????///////")






# #### tm lev density and scatter sampled

# print("tm density and scatter......")


# below_point_2_df = lev_tm_df.loc[lev_tm_df['tm_min'] <= 0.2]
# between_point_2_4 = lev_tm_df.loc[(lev_tm_df['tm_min'] > 0.2) & (lev_tm_df['tm_min'] <= 0.4) ]
# between_point_4_6 = lev_tm_df.loc[(lev_tm_df['tm_min'] > 0.4) & (lev_tm_df['tm_min'] <= 0.6) ]
# between_point_6_8 = lev_tm_df.loc[(lev_tm_df['tm_min'] > 0.6) & (lev_tm_df['tm_min'] <= 0.8) ]
# between_point_8_1 = lev_tm_df.loc[(lev_tm_df['tm_min'] > 0.8) & (lev_tm_df['tm_min'] <= 1) ]

# #between_point_6_8.to_csv("tm_between_point_6_8.csv")

# print(len(below_point_2_df),len(between_point_2_4),len(between_point_4_6),len(between_point_6_8),len(between_point_8_1))
# ##13380790 27122803 232548 4575 7162
# ##4019320 8227413 65897 2259 3277

# below_point_2_df = below_point_2_df.sample(2250)
# between_point_2_4 = between_point_2_4.sample(2250)
# between_point_4_6 = between_point_4_6.sample(2250)
# between_point_6_8 = between_point_6_8.sample(2000)
# between_point_8_1 = between_point_8_1.sample(2250)


# tm_all_scale = pd.concat([below_point_2_df,between_point_2_4,between_point_4_6,between_point_6_8,between_point_8_1], axis=0)

# # lev_dist = tm_all_scale['lev_distance'].to_list()
# # seq_dist = tm_all_scale['seq_idnt'].to_list()
# # tm_score = tm_all_scale["tm_min"].to_list()

# seq_dist_SS = tm_all_scale['lev%_SS'].to_list()
# seq_dist_AA = tm_all_scale['lev%_AA'].to_list()
# tm_score_min = tm_all_scale["tm_min"].to_list()
# tm_score_max = tm_all_scale["tm_max"].to_list()




# # below_50 = 0
# # between_50_70 = 0
# # above_70 = 0

# # good_tm = 0
# # bad_tm = 0


# # for i,value in enumerate(lev_dist):

# #     if tm_score[i] < 0.5:
# #         bad_tm = bad_tm + 1
# #     else:
# #         good_tm = good_tm + 1

# #     if lev_dist[i] < 50:
        
# #         below_50 = below_50 + 1
        
# #     elif lev_dist[i] > 50 and lev_dist[i] < 100:
        
# #         between_50_70 = between_50_70 + 1
    
# #     else:
        
# #         above_70 = above_70 + 1
        
# # print("lev below 50",below_50,"lev in between 50 to 100 ",between_50_70,"lev in between above 100 ", above_70)
# # print("tm below 0.5",bad_tm,"tm above 0.5",good_tm) 

# # print("tm density and scatter......done")

# seq_dist_SS = np.nan_to_num(seq_dist_SS)
# tm_score_min = np.nan_to_num(tm_score_min)
# seq_dist_AA = np.nan_to_num(seq_dist_AA)
# tm_score_max = np.nan_to_num(tm_score_max)


# # corr, p_v = pearsonr(lev_dist,tm_score)
# # print('Pearsons correlation Lev vs TM: %.3f' % corr, p_v)


# # temp_dict = {'tm_score':tm_score, 'lev_distance':lev_dist}

# # temp_df = pd.DataFrame(temp_dict)


# # print("kde plots lev_vs_TM")



# # g=sns.kdeplot(data=temp_df, x="tm_score", y="lev_distance",fill=True)
# # plt.savefig('tm_seq_idnt_den_scale_0-1_10k_original_SS_test_set_2.png', format="png")
# # plt.close()





# corr, p_v = pearsonr(seq_dist_SS,tm_score_min)
# print('Pearsons correlation TM-min vs lev_SS 0-1: %.3f' % corr, p_v)


# temp_dict = {'tm_score_min':tm_score_min, 'seq_ident_SS':seq_dist_SS}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_TM")



# g=sns.kdeplot(data=temp_df, x="tm_score_min", y="seq_ident_SS",fill=True)
# plt.savefig('tm_min_seq_ident_den_scale_0-1_10k_original_SS_test_set_2.png', format="png")
# plt.close()


# print("kde done lev_vs_TM")  

# ########
# corr, p_v = pearsonr(seq_dist_AA,tm_score_min)
# print('Pearsons correlation TM-min vs lev_AA 0-1: %.3f' % corr, p_v)


# temp_dict = {'tm_score_min':tm_score_min, 'seq_ident_AA':seq_dist_AA}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_TM min AA")



# g=sns.kdeplot(data=temp_df, x="tm_score_min", y="seq_ident_AA",fill=True)
# plt.savefig('tm_min_seq_ident_den_scale_0-1_10k_original_AA_test_set_2.png', format="png")
# plt.close()


# print("kde done lev_vs_TM min AA") 

# ########
# corr, p_v = pearsonr(seq_dist_SS,tm_score_max)
# print('Pearsons correlation TM-max vs lev_SS 0-1: %.3f' % corr, p_v)


# temp_dict = {'tm_score_max':tm_score_max, 'seq_ident_SS':seq_dist_SS}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_TM max SS")



# g=sns.kdeplot(data=temp_df, x="tm_score_max", y="seq_ident_SS",fill=True)
# plt.savefig('tm_max_seq_ident_den_scale_0-1_10k_original_SS_test_set_2.png', format="png")
# plt.close()


# print("kde done lev_vs_TM max SS") 

# ############

# corr, p_v = pearsonr(seq_dist_AA,tm_score_max)
# print('Pearsons correlation TM-max vs lev_AA 0-1: %.3f' % corr, p_v)


# temp_dict = {'tm_score_max':tm_score_max, 'seq_ident_AA':seq_dist_AA}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_TM max SS")



# g=sns.kdeplot(data=temp_df, x="tm_score_max", y="seq_ident_AA",fill=True)
# plt.savefig('tm_max_seq_ident_den_scale_0-1_10k_original_AA_test_set_2.png', format="png")
# plt.close()


# print("kde done lev_vs_TM max SAA") 


















