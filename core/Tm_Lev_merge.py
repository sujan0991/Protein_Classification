import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import gaussian_kde
from matplotlib.gridspec import GridSpec

# tm_df = pd.read_csv('merged.csv')
# lev_df = pd.read_csv('lev_distance_comparison_500.csv')



# tm_df.columns = ['domain1', 'domain2', 'length', 'rmsd', 'tm_query', 'tm_subject', 'cigar']



# new_df = pd.merge(lev_df, tm_df, how="outer", on=['domain1','domain2'])


# #print("len.........",len(new_df),new_df.head())

# new_df.to_csv('tm_lev.csv',index_label='index')#


#####################

tm_lev_df = pd.read_csv('tm_lev.csv')

lev_dis = tm_lev_df['lev_distance'].to_list()
tm_q = tm_lev_df['tm_query'].to_list()
tm_s = tm_lev_df['tm_subject'].to_list()

tm_mean_list = []

lev_dis_inverse_list = []

below_50 = 0
between_50_70 = 0
above_70 = 0

good_tm = 0
bad_tm = 0


# length = pd.Series(lev_dis) # convert to Series, so that it can be plotted
# sns.distplot(length, kde = False, 
# kde_kws = {'linewidth': 3},
# label = "distribution")

# plt.xlabel('distance', fontsize=14)
# plt.ylabel('Count', fontsize=14)
# #plt.savefig("500_lev_distribution.svg", format="svg")
# plt.show()




for i,value in enumerate(tm_q):
    
    
    
    tm_mean = (value + tm_s[i])/2
    
    if tm_mean < 0.5:
        bad_tm = bad_tm + 1
    else:
        good_tm = good_tm + 1
        
            
    
    tm_mean_list.append(tm_mean)
    
    if lev_dis[i] < 50:
        
        below_50 = below_50 + 1
        
    elif lev_dis[i] > 50 and lev_dis[i] < 70:
        
        between_50_70 = between_50_70 + 1
    
    else:
        
        above_70 = above_70 + 1
            
    
    
    
    
    
    
    if lev_dis[i] != 0:
        
        temp_lev = 1/lev_dis[i]
        lev_dis_inverse_list.append(temp_lev)
        
    else:
       
        lev_dis_inverse_list.append(lev_dis[i])
        
           
    
# length = pd.Series(tm_mean_list) # convert to Series, so that it can be plotted
# sns.distplot(length, kde = False, 
# kde_kws = {'linewidth': 3},
# label = "distribution")

# plt.xlabel('TM score', fontsize=14)
# plt.ylabel('Count', fontsize=14)
# #plt.savefig("500_TM_distribution.svg", format="svg")    
# plt.show()





# temp_dict = {'lev_distance':lev_dis, 'tm_score':tm_mean_list}

# temp_df = pd.DataFrame(temp_dict)

#print(temp_df.head())


# print(temp_df.corr(method ='pearson'))



#sns.scatterplot(data = temp_df, x = "lev_distance", y = "tm_score")

###########

# fig = plt.figure()
# gs = GridSpec(4, 4)

# ax_scatter = fig.add_subplot(gs[1:4, 0:3])
# ax_hist_x = fig.add_subplot(gs[0,0:3])
# ax_hist_y = fig.add_subplot(gs[1:4, 3])

# ax_scatter.scatter(temp_df['lev_distance'], temp_df['tm_score'])

# ax_hist_x.hist(temp_df['lev_distance'])
# ax_hist_y.hist(temp_df['tm_score'], orientation = 'horizontal')

# # plt.xlabel('Levenshtein Distance (Sequential)', fontsize=10)
# # plt.ylabel('TM score (Structural)', fontsize=10)
# plt.show()


#############

# lev_dis = np.nan_to_num(lev_dis)
# tm_mean_list = np.nan_to_num(tm_mean_list)

temp_dict = {'lev_distance':lev_dis, 'tm_score':tm_mean_list}

temp_df = pd.DataFrame(temp_dict)

# x = temp_df['lev_distance']
# y = temp_df['tm_score']

# # Calculate the point density
# xy = np.vstack([x,y])
# z = gaussian_kde(xy)(xy)

# fig, ax = plt.subplots()
# ax.scatter(x, y, c=z, s=100)
# plt.show()





# fig, ax = plt.subplots(figsize=(6, 6))
# fig, ax = plt.subplots(figsize=(6, 6))
# sns.scatterplot(
#     data=temp_df,
#     x="lev_distance",
#     y="tm_score",
#     color= "k",
#     ax=ax,
# )
# sns.kdeplot(
#     data=temp_df,
#     x="lev_distance",
#     y="tm_score",
#     levels=5,
#     fill=True,
#     alpha=0.6,
#     cut=2,
#     ax=ax,
# )

# plt.show()

# print("kde plots")
# g=sns.kdeplot(data=temp_df, x="lev_distance", y="tm_score",fill=True)
# plt.show()
# #plt.savefig(os.path.join(outdir, "blast_kde.png"), transparent=True, bbox_inches="tight", dpi=300)
# #plt.close()



# exit()

############

# g2 = sns.lineplot(data = temp_df, x = "lev_distance", y = "tm_score",ci=None)
# sns.despine()
# g2.set(xlabel=None)
# g2.set(ylabel=None)
# plt.xlabel('TM score (Structural)', fontsize=14)
# plt.ylabel('Levenshtein Distance (Sequential)', fontsize=14)
# plt.show()


fig = plt.figure()
gs = GridSpec(4, 4)

ax_scatter = fig.add_subplot(gs[1:4, 0:3])
ax_hist_x = fig.add_subplot(gs[0,0:3])
ax_hist_y = fig.add_subplot(gs[1:4, 3])

ax_scatter.lines(temp_df['lev_distance'], temp_df['tm_score'])


ax_hist_x.hist(temp_df['lev_distance'])
ax_hist_y.hist(temp_df['tm_score'], orientation = 'horizontal')

# plt.xlabel('Levenshtein Distance (Sequential)', fontsize=10)
# plt.ylabel('TM score (Structural)', fontsize=10)
plt.show()

exit()
    

print(below_50,between_50_70,above_70)
print("tm",good_tm,bad_tm)
#print(len(tm_mean_list))      



print(np.isnan(lev_dis).any(),np.isinf(lev_dis).any())

lev_dis = np.nan_to_num(lev_dis)

lev_dis_inverse_list = np.nan_to_num(lev_dis_inverse_list)

corr, p_v = pearsonr(tm_mean_list,lev_dis)
print('Pearsons correlation TM vs Lev: %.3f' % corr, p_v)

# corr, _ = spearmanr(tm_mean_list,lev_dis)
# print('Spearmans correlation TM vs Lev: %.3f' % corr)

corr, p_v = pearsonr(tm_mean_list,lev_dis_inverse_list)
print('Pearsons correlation.......TM vs 1/Lev: %.3f' % corr, p_v)

# corr, _ = spearmanr(tm_mean_list,lev_dis_inverse_list)
# print('Spearmans correlation:.........TM vs 1/Lev %.3f' % corr)


temp_dict = {'1/lev_distance':lev_dis_inverse_list, 'tm_score':tm_mean_list}

temp_df = pd.DataFrame(temp_dict)

#print(temp_df.head())


# print(temp_df.corr(method ='pearson'))



#sns.scatterplot(data = temp_df, x = "lev_distance", y = "tm_score")

# sns.lineplot(data = temp_df, x = "1/lev_distance", y = "tm_score",ci=None).set(title='TM_Score vs 1/Lev_Distance')

# plt.show()
