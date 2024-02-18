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

# lev_df = pd.read_csv('distance_comparison_test_set(0.2)_10k_superfamily_domain_without_AA.csv')
# tm_df = pd.read_csv('tm_merge_10k_test_set_2.csv')

# # print("",lev_df.columns.tolist())

# lev_df.columns = ['index', 'domain1', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'lev_distance', 'domain2', 'Class2', 'Arch2', 'Topol2', 'Homol2']
# tm_df.columns = ['domain1', 'domain2', 'length', 'rmsd', 'tm_query', 'tm_subject', 'cigar']

# print("lev_df len",len(lev_df),len(tm_df))

# tm_not_same_domain = tm_df[(tm_df['domain1'] != tm_df['domain2'])]

# lev_df['domains'] = ['-'.join(sorted(combine)) for combine in zip(lev_df['domain1'], lev_df['domain2'])]

# tm_df['domains'] = ['-'.join(sorted(combine)) for combine in zip(tm_df['domain1'], tm_df['domain2'])]


# lev_tm_df = lev_df.merge(tm_df,left_on='domains',right_on='domains')

# print("///////",lev_tm_df.columns.tolist())

# lev_tm_df = lev_tm_df[lev_tm_df['lev_distance'].notna()]
# lev_tm_df = lev_tm_df[lev_tm_df['tm_query'].notna()]

# print("lev_tm_df len",len(lev_tm_df))


# lev_tm_df.to_csv('lev_tm_10k_without_AA_2_df.csv')


# exit()

##############






# lev_tm_df = pd.read_csv('lev_tm_10k_without_AA_with_seq_df_2.csv')

# lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)


#print(len(lev_tm_df),lev_tm_df.columns.to_list())


# ##### lev vs lev of original AA
# lev_tm_original_AA_df = pd.read_csv('lev_tm_AA_10k_df.csv')

# lev_dist = lev_tm_df['lev_distance'].to_list()
# lev_dist_AA = lev_tm_original_AA_df["lev_distance_AA"].to_list()

# lev_dist = np.nan_to_num(lev_dist)
# lev_dist_AA = np.nan_to_num(lev_dist_AA)
# lev_dist = lev_dist[:len(lev_dist_AA)]


# corr, p_v = pearsonr(lev_dist,lev_dist_AA)
# print('Pearsons correlation Lev vs lev_dist_AA: %.3f' % corr, p_v)


# temp_dict = {'lev_dist_AA':lev_dist_AA, 'lev_distance':lev_dist}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_TM")



# g=sns.lineplot(data=temp_df, x="lev_dist_AA", y="lev_distance")
# plt.savefig('lev_without_AA_lev_original_AA_line_10k.png', format="png")
# plt.close()


# exit()



#### seq identity plot ## see in density plot section

lev_tm_df['seq_1_len'] = lev_tm_df['sequence_1'].apply(len)
 
lev_tm_df['seq_2_len'] = lev_tm_df['sequence_2'].apply(len)
 
lev_tm_df['max_seq_len'] = lev_tm_df[["seq_1_len", "seq_2_len"]].max(axis=1)
 
lev_tm_df['seq_idnt'] = lev_tm_df['lev_distance']/ lev_tm_df['max_seq_len']

# #

# ##### lev, tm distribution plot

# # print("lev, tm distribution plot")




# # below_50 = 0
# # between_50_70 = 0
# # above_70 = 0

# # good_tm = 0
# # bad_tm = 0

# # lev_dist = lev_tm_df['lev_distance'].to_list()
# # tm_score = lev_tm_df["tm_min"].to_list()

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

# # print("lev, tm distribution plot done")

# # sns.displot(lev_tm_df["lev_distance"])
# # plt.savefig('lev_distribution_10k_without_AA_test_set_2.png', format="png")
# # plt.close()

# # sns.displot(lev_tm_df["tm_min"])
# # plt.savefig('tm_distribution_10k_without_AA_test_set_2.png', format="png")
# # plt.close()


# # exit()






# ########tm lev density


# # lev_dist = lev_tm_df['lev_distance'].to_list()
# # tm_score = lev_tm_df["tm_min"].to_list()

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

# # lev_dist = np.nan_to_num(lev_dist)
# # tm_score = np.nan_to_num(tm_score)


# # corr, p_v = pearsonr(lev_dist,tm_score)
# # print('Pearsons correlation Lev vs tm: %.3f' % corr, p_v)


# # temp_dict = {'tm_score':tm_score, 'lev_distance':lev_dist}

# # temp_df = pd.DataFrame(temp_dict)


# # print("kde plots lev_vs_TM")



# # g=sns.kdeplot(data=temp_df, x="tm_score", y="lev_distance",fill=True)
# # plt.savefig('tm_lev_den_10k_without_AA.png', format="png")
# # plt.close()








# #### tm lev density and scatter sampled

# print("tm density and scatter......")



# lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)



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
# between_point_6_8 = between_point_6_8.sample(2250)
# between_point_8_1 = between_point_8_1.sample(2250)


# tm_all_scale = pd.concat([below_point_2_df,between_point_2_4,between_point_4_6,between_point_6_8,between_point_8_1], axis=0)

# lev_dist = tm_all_scale['lev_distance'].to_list()
# seq_dist = tm_all_scale['seq_idnt'].to_list()
# tm_score = tm_all_scale["tm_min"].to_list()

# below_50 = 0
# between_50_70 = 0
# above_70 = 0

# good_tm = 0
# bad_tm = 0


# for i,value in enumerate(lev_dist):

#     if tm_score[i] < 0.5:
#         bad_tm = bad_tm + 1
#     else:
#         good_tm = good_tm + 1

#     if lev_dist[i] < 50:
        
#         below_50 = below_50 + 1
        
#     elif lev_dist[i] > 50 and lev_dist[i] < 100:
        
#         between_50_70 = between_50_70 + 1
    
#     else:
        
#         above_70 = above_70 + 1
        
# print("lev below 50",below_50,"lev in between 50 to 100 ",between_50_70,"lev in between above 100 ", above_70)
# print("tm below 0.5",bad_tm,"tm above 0.5",good_tm) 

# print("tm density and scatter......done")

# lev_dist = np.nan_to_num(lev_dist)
# tm_score = np.nan_to_num(tm_score)


# corr, p_v = pearsonr(lev_dist,tm_score)
# print('Pearsons correlation Lev vs TM: %.3f' % corr, p_v)


# temp_dict = {'tm_score':tm_score, 'lev_distance':lev_dist}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_TM")



# g=sns.kdeplot(data=temp_df, x="tm_score", y="lev_distance",fill=True)
# plt.savefig('tm_lev_distance_den_scale_0-1_10k_without_AA_test_set_2.png', format="png")
# plt.close()



# corr, p_v = pearsonr(seq_dist,tm_score)
# print('Pearsons correlation TM vs seq. Idnt: %.3f' % corr, p_v)


# temp_dict = {'tm_score':tm_score, 'seq_ident':seq_dist}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_TM")



# g=sns.kdeplot(data=temp_df, x="tm_score", y="seq_ident",fill=True)
# plt.savefig('tm_seq_ident_den_scale_0-1_10k_without_AA_test_set_2.png', format="png")
# plt.close()


# print("kde done lev_vs_TM")  

# sns.scatterplot(data= temp_dict, x="tm_score", y="lev_distance")
# plt.savefig('tm_lev_scater_scale_0-1_10k_without_AA.png', format="png")
# plt.close()



# exit()

# ########### lev superfamily

# same_homo = lev_tm_df.query('Homol1 == Homol2') ['lev_distance']
# diff_homo = lev_tm_df.query('Homol1 != Homol2') ['lev_distance']

# same_homo_df = pd.DataFrame({'lev_distance':same_homo.values, 'Catagory':'Same Superfamily'})
# diff_homo_df = pd.DataFrame({'lev_distance':diff_homo.values, 'Catagory':'Different Superfamily'})
 
# same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

# g2 = sns.violinplot( y=same_homo_not_same_homo_df["lev_distance"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# #g2.set(ylabel=None)
# plt.savefig('lev_super_10k_without_AA_test_set_2_v.png', format="png")
# plt.close()



# # ### tm superfamily plot

# # # print("tm superfamily plot......")
# # # lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)

# # # same_homo = lev_tm_df.query('Homol1 == Homol2') ['tm_min']
# # # not_same_homo = lev_tm_df.query('Homol1 != Homol2') ['tm_min']


# # # below = same_homo[same_homo <= 0.5]

# # # above = same_homo[same_homo > 0.5]

# # # print("same same_homo len",len(same_homo),"below 0.5",len(below),"above 0.5",len(above))



# # # below_n = not_same_homo[not_same_homo <= 0.5]

# # # above_n = not_same_homo[not_same_homo > 0.5]

# # # print("not same homo len",len(not_same_homo),"below 0.5",len(below_n),"above 0.5",len(above_n))


# # # same_homo__df = pd.DataFrame({'TM Score':same_homo.values, 'Catagory':'Same Superfamily'})
# # # diff_homo_df = pd.DataFrame({'TM Score':not_same_homo.values, 'Catagory':'Different Superfamily'})
 
# # # same_homo_not_same_homo_above_point_5_df = same_homo__df.append(diff_homo_df, ignore_index=True)

# # # g2 = sns.violinplot( y=same_homo_not_same_homo_above_point_5_df["TM Score"], x=same_homo_not_same_homo_above_point_5_df["Catagory"],bw=.05 )
# # # sns.despine()
# # # g2.set(xlabel=None)
# # # #g2.set(ylabel=None)
# # # plt.savefig('lev_super_below_tm_5_10k_c-12-26_v_resi_100-200.png', format="png")
# # # plt.close()
# # # exit() 





# # ## tm above and below 0.5 for lev

# print("lev........")


# below_point_5= lev_tm_df[lev_tm_df.tm_min<0.5]
# abobe_point_5= lev_tm_df[lev_tm_df.tm_min>=0.5]
    


# below_point_5_same_homo = below_point_5.query('Homol1 == Homol2') ['lev_distance']
# below_point_5_diff_homo = below_point_5.query('Homol1 != Homol2') ['lev_distance']

# above_point_5_same_homo = abobe_point_5.query('Homol1 == Homol2') ['lev_distance']
# above_point_5_diff_homo = abobe_point_5.query('Homol1 != Homol2') ['lev_distance']
    

# ## below plot
# same_homo_below_point_5_df = pd.DataFrame({'lev_distance':below_point_5_same_homo.values, 'Catagory':'Same Superfamily'})
# diff_homo_below_point_5_df = pd.DataFrame({'lev_distance':below_point_5_diff_homo.values, 'Catagory':'Different Superfamily'})
 
# same_homo_not_same_homo_below_point_5_df = same_homo_below_point_5_df.append(diff_homo_below_point_5_df, ignore_index=True)

# g2 = sns.violinplot( y=same_homo_not_same_homo_below_point_5_df["lev_distance"], x=same_homo_not_same_homo_below_point_5_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# #g2.set(ylabel=None)
# plt.savefig('lev_super_below_tm_5_10k_without_AA_test_set_2_v.png', format="png")
# plt.close()
# ## above plot
# same_homo_above_point_5_df = pd.DataFrame({'lev_distance':above_point_5_same_homo.values, 'Catagory':'Same Superfamily'})
# diff_homo_above_point_5_df = pd.DataFrame({'lev_distance':above_point_5_diff_homo.values, 'Catagory':'Different Superfamily'})
 
# same_homo_not_same_homo_above_point_5_df = same_homo_above_point_5_df.append(diff_homo_above_point_5_df, ignore_index=True)

# g2 = sns.violinplot( y=same_homo_not_same_homo_above_point_5_df["lev_distance"], x=same_homo_not_same_homo_above_point_5_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# #g2.set(ylabel=None)
# plt.savefig('lev_super_above_tm_5_10k_without_AA_test_set_2_v.png', format="png")
# plt.close()

# print("tm superfamily plot......done")


# ##################






# ####### violin plot all

# print("violin plot all")

# lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)

# lev_tm_df["tm_min"] = lev_tm_df["tm_min"].round(decimals = 1)


# lev_dist = lev_tm_df['lev_distance'].to_list()
# tm_score = lev_tm_df["tm_min"].to_list()


# tm_0 = 0
# tm_1 =0
# tm_2 = 0
# tm_3 = 0
# tm_4 =0
# tm_5 = 0
# tm_6 = 0
# tm_7 =0
# tm_8 = 0
# tm_9 = 0
# tm_10 =0


# for i,value in enumerate(tm_score):

#     if tm_score[i] == 0.0:
#         tm_0 = tm_0 + 1
#     elif tm_score[i] == 0.1:
#         tm_1 = tm_1 + 1 
#     elif tm_score[i] == 0.2:
#         tm_2 = tm_2 + 1
        
#     elif tm_score[i] == 0.3:
#         tm_3 = tm_3 + 1
        
#     elif tm_score[i] == 0.4:
#         tm_4 = tm_4 + 1
#     elif tm_score[i] == 0.5:
#         tm_5 = tm_5 + 1
#     elif tm_score[i] == 0.6:
#         tm_6 = tm_6 + 1
#     elif tm_score[i] == 0.7:
#         tm_7 = tm_7 + 1
#     elif tm_score[i] == 0.8:
#         tm_8 = tm_8 + 1  
#     elif tm_score[i] == 0.9:
#         tm_9 = tm_9 + 1 
#     else :                               
#         tm_10 = tm_10 + 1
        
# print("tm 0-1",tm_0,tm_1,tm_2,tm_3,tm_4,tm_5,tm_6,tm_7,tm_8,tm_9,tm_10) 

      
# below_50 = 0
# between_50_70 = 0
# above_70 = 0

# good_tm = 0
# bad_tm = 0


# for i,value in enumerate(lev_dist):

#     if tm_score[i] < 0.5:
#         bad_tm = bad_tm + 1
#     else:
#         good_tm = good_tm + 1

#     if lev_dist[i] < 50:
        
#         below_50 = below_50 + 1
        
#     elif lev_dist[i] > 50 and lev_dist[i] < 100:
        
#         between_50_70 = between_50_70 + 1
    
#     else:
        
#         above_70 = above_70 + 1
        
# print("lev below 50",below_50,"lev in between 50 to 100 ",between_50_70,"lev in between above 100 ", above_70)
# print("tm below 0.5",bad_tm,"tm above 0.5",good_tm) 


# temp_dict = {'tm_score':tm_score, 'lev_distance':lev_dist}

# temp_df = pd.DataFrame(temp_dict)

# print("temp_df........len",len(temp_df))

# sns.violinplot(data=temp_df, x='tm_score', y='lev_distance')
# plt.savefig('tm_lev_violin_all_without_AA_test_set_2.png', format="png")
# plt.close()


# print("violin plot done")







####### when tm score 0.0

lev_tm_df["tm_min"] = lev_tm_df[["tm_query", "tm_subject"]].min(axis=1)

lev_tm_df["tm_min"] = lev_tm_df["tm_min"].round(decimals = 1)

tm_0 = lev_tm_df[lev_tm_df.tm_min == 0.0]


same_homo = tm_0.query('Homol1 == Homol2') ['lev_distance']
diff_homo = tm_0.query('Homol1 != Homol2') ['lev_distance']


print("same superfaamily,diff superfamily",len(same_homo),len(diff_homo))

# same_homo_df = pd.DataFrame({'lev_distance':same_homo.values, 'Catagory':'Same Superfamily'})
# diff_homo_df = pd.DataFrame({'lev_distance':diff_homo.values, 'Catagory':'Different Superfamily'})
 
# same_homo_not_same_homo_df = same_homo_df.append(diff_homo_df, ignore_index=True)

# g2 = sns.violinplot( y=same_homo_not_same_homo_df["lev_distance"], x=same_homo_not_same_homo_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# #g2.set(ylabel=None)
# plt.savefig('tm_0_10k_without_AA_test_set_2_v.png', format="png")
# plt.close()



##############

# group1 = lev_tm_df.query('tm_min < 0.5') ['lev_distance']
# group2 = lev_tm_df.query('tm_min >= 0.5') ['lev_distance']

# group1_df = pd.DataFrame({'lev_distance':group1.values, 'Catagory':'tm_min < 0.4'})
# group2_df = pd.DataFrame({'lev_distance':group2.values, 'Catagory':'tm_min >= 0.4'})

# group1_group2_df = group1_df.append(group2_df, ignore_index=True)

# g2 = sns.violinplot( y=group1_group2_df["lev_distance"], x=group1_group2_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# plt.savefig('tm_lev_violin_below_above_5_without_AA_test_set_2.png', format="png")
# plt.close() 




########## 

# print("???????????")
# group1 = lev_tm_df.query('tm_min >= 0.5 & Homol1 == Homol2') ['lev_distance']
# group2 = lev_tm_df.query('tm_min < 0.5 & Homol1 == Homol2') ['lev_distance']
# group3 = lev_tm_df.query('tm_min >= 0.5 & Homol1 != Homol2') ['lev_distance']
# group4 = lev_tm_df.query('tm_min < 0.5 & Homol1 != Homol2') ['lev_distance']


# group1_df = pd.DataFrame({'lev_distance':group1.values, 'Catagory':'tm_min >= 0.5 & same S.F'})
# group2_df = pd.DataFrame({'lev_distance':group2.values, 'Catagory':'tm_min < 0.5 & same S.F'})
# group3_df = pd.DataFrame({'lev_distance':group3.values, 'Catagory':'tm_min >= 0.5 & diff. S.F'})
# group4_df = pd.DataFrame({'lev_distance':group4.values, 'Catagory':'tm_min < 0.5 & diff. S.F'})

# print(len(group1_df),len(group2_df),len(group3_df),len(group4_df))

# group1_group2_df = pd.concat([group1_df,group2_df,group3_df,group4_df])


# g2 = sns.violinplot( y=group1_group2_df["lev_distance"], x=group1_group2_df["Catagory"],bw=.05 )
# sns.despine()
# g2.set(xlabel=None)
# plt.savefig('tm_lev_violin_below_above_5_same_diff_super_without_AA_test_set_2.png', format="png")
# plt.close() 

# print("???????????///////")









exit()

































# ## below 50
#below_50_df = pd.read_csv("lev_below_50_test_5k.csv")


#tm_not_same_domain = swap_columns(tm_not_same_domain, 'domain1', 'domain2')

tm_above_point_5_df= above_point_5_df.merge(lev_df, on=['domain1','domain2'])

print("merged tm_above_point_5_df len.........",len(tm_above_point_5_df))


tm_above_point_5_df = tm_above_point_5_df[tm_above_point_5_df['lev_distance'].notna()]


print("............",len(tm_above_point_5_df))


tm_above_point_5_df.to_csv('tm_priority_lev_50k_df.csv')




exit()



above_point_5_df = pd.read_csv("tm_min_below_point_5_df_test_5k.csv")

print("",lev_df.columns.tolist())

lev_df.columns = ['index', 'domain1', 'Class1', 'Arch1', 'Topol1', 'Homol1', 'lev_distance', 'domain2', 'Class2', 'Arch2', 'Topol2', 'Homol2']

#tm_not_same_domain = tm_df[(tm_df['domain1'] != tm_df['domain2'])]


# ## below 50
#below_50_df = pd.read_csv("lev_below_50_test_5k.csv")


#tm_not_same_domain = swap_columns(tm_not_same_domain, 'domain1', 'domain2')

tm_above_point_5_df= above_point_5_df.merge(lev_df, on=['domain1','domain2'])

print("merged tm_above_point_5_df len.........",len(tm_above_point_5_df))


tm_above_point_5_df = tm_above_point_5_df[tm_above_point_5_df['lev_distance'].notna()]


print("............",len(tm_above_point_5_df))

tm_above_point_5_df.to_csv('tm_min_below_point_5_lev_df.csv')



exit()



below_50 = 0
between_50_70 = 0
above_70 = 0

good_tm = 0
bad_tm = 0

lev_tm_df = pd.read_csv('tm_lev_10k.csv')

lev_tm_df = lev_tm_df.sample(frac=0.001).reset_index(drop=True)

print("3k len",len(lev_tm_df))

lev_dist = lev_tm_df['lev_distance'].to_list()

tm_q = lev_tm_df['tm_query'].to_list()
tm_s = lev_tm_df['tm_subject'].to_list()

tm_mean_list = []

for i,value in enumerate(tm_q):

    tm_mean = (value + tm_s[i])/2
    
    tm_mean_list.append(tm_mean)


    if tm_mean < 0.5:
        bad_tm = bad_tm + 1
    else:
        good_tm = good_tm + 1

    if lev_dist[i] < 50:
        
        below_50 = below_50 + 1
        
    elif lev_dist[i] > 50 and lev_dist[i] < 100:
        
        between_50_70 = between_50_70 + 1
    
    else:
        
        above_70 = above_70 + 1
        
        
print("tm_mean_list len",len(tm_mean_list),len(lev_dist))

print("lev below 50",below_50,"lev in between 50 to 100 ",between_50_70,"lev in between above 100 ", above_70)
print("tm below 0.5",bad_tm,"tm above 0.5",good_tm)   


exit()


tm_mean_list = np.nan_to_num(tm_mean_list)
lev_dist = np.nan_to_num(lev_dist)

corr, p_v = pearsonr(lev_dist,tm_mean_list)
print('Pearsons correlation TM vs lev: %.3f' % corr, p_v)




temp_dict = {'tm_score':tm_mean_list, 'Lev_distance':lev_dist}

temp_df = pd.DataFrame(temp_dict)



# # sns.lineplot(data = temp_df, x = "Lev_distance", y = "tm_score",ci=None)

# # plt.savefig('Lev_vs_TM_line_10k_3k_c-12-26.svg', format="svg")
# # #plt.show()



print("kde plots lev_vs_TM")

g=sns.kdeplot(data=temp_df, x="Lev_distance", y="tm_score",fill=True)
plt.savefig('lev_vs_TM_den_10k_3k_c-12-26.svg', format="svg")
#plt.show()
print("kde done lev_vs_TM")     
