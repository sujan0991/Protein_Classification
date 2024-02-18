import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from scipy.stats import pearsonr


######### plot #########

lev_tm_original_AA_df = pd.read_csv('lev_tm_AA_10k_df.csv')

print(lev_tm_original_AA_df.head())


#### tm lev density 

print("tm density ......")



# lev_tm_original_AA_df["tm_min"] = lev_tm_original_AA_df[["tm_query", "tm_subject"]].min(axis=1)



# below_point_2_df = lev_tm_original_AA_df.loc[lev_tm_original_AA_df['tm_min'] <= 0.2]
# between_point_2_4 = lev_tm_original_AA_df.loc[(lev_tm_original_AA_df['tm_min'] > 0.2) & (lev_tm_original_AA_df['tm_min'] <= 0.4) ]
# between_point_4_6 = lev_tm_original_AA_df.loc[(lev_tm_original_AA_df['tm_min'] > 0.4) & (lev_tm_original_AA_df['tm_min'] <= 0.6) ]
# between_point_6_8 = lev_tm_original_AA_df.loc[(lev_tm_original_AA_df['tm_min'] > 0.6) & (lev_tm_original_AA_df['tm_min'] <= 0.8) ]
# between_point_8_1 = lev_tm_original_AA_df.loc[(lev_tm_original_AA_df['tm_min'] > 0.8) & (lev_tm_original_AA_df['tm_min'] <= 1) ]

# between_point_6_8.to_csv("tm_between_point_6_8_original_AA.csv")

# print(len(below_point_2_df),len(between_point_2_4),len(between_point_4_6),len(between_point_6_8),len(between_point_8_1))
# ##13380790 27122803 232548 4575 7162


# below_point_2_df = below_point_2_df.sample(4500)
# between_point_2_4 = between_point_2_4.sample(4500)
# between_point_4_6 = between_point_4_6.sample(4500)
# between_point_6_8 = between_point_6_8.sample(4166)
# between_point_8_1 = between_point_8_1.sample(4500)


# tm_all_scale = pd.concat([below_point_2_df,between_point_2_4,between_point_4_6,between_point_6_8,between_point_8_1], axis=0)

# lev_dist = tm_all_scale['lev_distance_AA'].to_list()
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
# print('Pearsons correlation Lev vs tm original AA: %.3f' % corr, p_v)



# temp_dict = {'tm_score':tm_score, 'lev_distance':lev_dist}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_lev original AA")

# g=sns.kdeplot(data=temp_df, x="tm_score", y="lev_distance",fill=True)
# plt.savefig('tm_lev_den_scale_0-1_50k_original_AA.png', format="png")
# plt.close()



####### our method lev vs original AA lev  ########


# lev_dist_AA = lev_tm_original_AA_df['lev_distance_AA'].to_list()
# lev_dist = lev_tm_original_AA_df["lev_distance"].to_list()


# lev_dist = np.nan_to_num(lev_dist)
# lev_dist_AA = np.nan_to_num(lev_dist_AA)

# corr, p_v = pearsonr(lev_dist,lev_dist_AA)
# print('Pearsons correlation Lev vs lev original AA: %.3f' % corr, p_v)


# temp_dict = {'lev_distance_AA':lev_dist_AA, 'lev_distance':lev_dist}

# temp_df = pd.DataFrame(temp_dict)


# print("kde plots lev_vs_lev original AA done")

# # g=sns.kdeplot(data=temp_df, x="lev_distance_AA", y="lev_distance",fill=True)
# # g=sns.kdeplot(data=temp_df, x="lev_distance_AA", y="lev_distance")
# # plt.savefig('lev_den_scale_0-1_50k_original_AA.png', format="png")
# # plt.close()

# sns.scatterplot(data=temp_df, x="lev_distance_AA", y="lev_distance")
# plt.savefig('lev_vs_lev_scatter_original_AA_.png', format="png")
# plt.close()





tm_between_point_6_8_df = pd.read_csv("tm_between_point_6_8_original_AA.csv")

print(len(tm_between_point_6_8_df))

# group1 = tm_between_point_6_8_df[tm_between_point_6_8_df.lev_distance>130]
# group2 = tm_between_point_6_8_df[tm_between_point_6_8_df.lev_distance<130]

group1 = tm_between_point_6_8_df.query('lev_distance > 130') ['lev_distance_AA']
group2 = tm_between_point_6_8_df.query('lev_distance < 130') ['lev_distance_AA']

group1_df = pd.DataFrame({'lev_distance_AA':group1.values, 'Catagory':'Lev > 130'})
group2_df = pd.DataFrame({'lev_distance_AA':group2.values, 'Catagory':'Lev < 130'})
 
group1_group2_df = group1_df.append(group2_df, ignore_index=True)

# print(group1_group2_df.head())


g2 = sns.violinplot( y=group1_group2_df["lev_distance_AA"], x=group1_group2_df["Catagory"],bw=.05 )
sns.despine()
g2.set(xlabel=None)
#g2.set(ylabel=None)
# plt.savefig('lev_super_below_tm_5_10k_c-12-26_v.png', format="png")
plt.show()