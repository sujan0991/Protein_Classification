## This script is to do statistical analysis of the features
# import pymol
# import pymol.cmd as cmd

# import psico.editing

import pandas as pd
import sys
# from pymol import stored
import multiprocessing as mp
import os.path
import pickle
import statistics
from statistics import mode
from collections import Counter
import matplotlib.pyplot as plt
import random
import seaborn as sns
import numpy as np
from scipy.stats import genextreme


# Return the length of feature vector. As all features in one vector have same length,
# returning one's length is ok. Here returning len(temp_dic['resns']

def length_calculation(values):
    # mp.current_process()
    temp_dic = {}
    temp_dic = values[1]
    # print("length of list",  len(temp_dic['resns']),len(temp_dic['ss_struc']),len(temp_dic['distances']))
    # if len(temp_dic["resns"]) > 53 and len(temp_dic["resns"]) < 249:
    #     print("len less then 249,> 54",values[0])   
    if  len(temp_dic["resns"]) < 3:
            print('len less than 3')
    return len(temp_dic["resns"])


def all_AA_list(values):
    temp_dic = {}
    temp_dic = values
    return temp_dic["resns"]


def all_ss_list(values):
    temp_dic = {}
    temp_dic = values
    return temp_dic["ss_struc"]


def all_distance_list(values):
    temp_dic = {}
    temp_dic = values
    return temp_dic["distances"]


aa_dict = {}

def aa_number_plot(aa_list):

    # get every aa count from aa_list
    aa_count = Counter(aa_list)

    print("aa_count.....LEU", aa_count["LEU"])

    for k, v in aa_count.items():
      if k == 'ALA':
          aa_dict[k] = v
      elif k == 'ARG':
          aa_dict[k] = v
      elif k == 'ASN':
          aa_dict[k] = v
      elif k == 'ASP':
          aa_dict[k] = v
      elif k == 'CYS':
          aa_dict[k] = v
      elif k == 'GLN':
          aa_dict[k] = v
      elif k == 'GLU':
          aa_dict[k] = v
      elif k == 'GLY':
          aa_dict[k] = v
      elif k == 'HIS':
          aa_dict[k] = v
      elif k == 'ILE':
          aa_dict[k] = v
      elif k == 'LEU':
          aa_dict[k] = v
      elif k == 'LYS':
          aa_dict[k] = v
      elif k == 'MET':
          aa_dict[k] = v
      elif k == 'PHE':
          aa_dict[k] = v
      elif k == 'PRO':
          aa_dict[k] = v
      elif k == 'SER':
          aa_dict[k] = v
      elif k == 'THR':
          aa_dict[k] = v
      elif k == 'TRP':
          aa_dict[k] = v
      elif k == 'TYR':
          aa_dict[k] = v
      elif k == 'VAL':
          aa_dict[k] = v
      
        

    plt.bar(list(aa_dict.keys()), aa_dict.values())
    plt.savefig('aa_count.png')
    #plt.show()


def ss_number_plot(ss_list):
    # get every ss count from ss_list
    ss_count = Counter(ss_list)
    print("ss_count.....H", ss_count["H"])
    plt.bar(ss_count.keys(), ss_count.values())
    plt.savefig('ss_count.png')
    #plt.show()




# conver list of lists to single list
def flatten(t):
    return [item for sublist in t for item in sublist]


if __name__ == "__main__":
    with open("/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/pickle_folder/cath_all_proteins_dictionary.pickle", "rb") as handle:
        all_proteins_dictionary = pickle.load(handle)
    print("all_proteins_dictionary len",len(all_proteins_dictionary))
    # print(dict(list(all_proteins_dictionary.items())[0:2]))
    # # run length_calculation in parallel for all PDBs using multiprocessing and store the results in length_of_dict_value
    with mp.Pool(mp.cpu_count()) as pool:
        length_of_dict_value = pool.map(
            length_calculation, all_proteins_dictionary.items()
        )
    newlist = list(filter(None, length_of_dict_value))
    print("length_of_final data",len(length_of_dict_value),len(newlist))
    df = pd.DataFrame(newlist, columns=['length'])
    # plot the distribution
sns.displot(df, x="length", kde=True)
plt.xlabel('Number of Residues')
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_all_cath.png', format="png")
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/final_plots/length_distribution_all_cath.svg', format="svg")
plt.close()
    # plt.savefig('data_dis.png')
    #plt.show()
    
    #####




test = dict(list(all_proteins_dictionary.items())[0:2])

for key, value in all_proteins_dictionary.items():
    temp_dic = {}
    temp_dic = value
    if  len(temp_dic["resns"]) < 3:
            print('len less than 3',key,temp_dic["resns"])
    
    
    # calculate the mean.,median, min, max  length of length_of_dict_value
    # mean = statistics.mean(newlist)
    # median = statistics.median(newlist)

    # print(
    #     "vector per protein....mean.,median. min. max....",
    #     mean,
    #     median,
        
    # )

    # # # run all_AA_list in parallel for all PDBS using multiprocessing and store the results in AA_list

    with mp.Pool(mp.cpu_count()) as pool:

        AA_list = pool.map(all_AA_list, all_proteins_dictionary.values())

    with mp.Pool(mp.cpu_count()) as pool:

        ss_list = pool.map(all_ss_list, all_proteins_dictionary.values())

    
    # convert list of lists to single list
    aa = flatten(AA_list)
    ss = flatten(ss_list)
    
    # Plot
    aa_number_plot(aa)
    ss_number_plot(ss)
   








# for now pick 0.1% of the datom randomly with seed 1111
random.seed(1111)
p = 1 # 0.1% of the lines

# path on workstation
path = '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_SS_BLAST_s20_results.csv'
data = pd.read_csv(path, header=0, skiprows=lambda i: i>0 and random.random() > p, usecols = ['cath_superFamily', 'ss_score'], dtype={'cath_superFamily': bool, 'ss_score': 'float64'}, engine='c', memory_map=True)


data_same = data[data.cath_superFamily]
data_diff = data[~data.cath_superFamily]

SS_scores = np.linspace(0.0,1.0,21)
delta = SS_scores[1]-SS_scores[0]
N_SS = np.zeros_like(SS_scores)
Nbar_SS = np.zeros_like(SS_scores)

# count up all pairs that fall in the interval bins
for i in range(len(SS_scores)):
    N_SS[i] = data_same[((SS_scores[i]-data_same.ss_score)>=0) & ((SS_scores[i]-data_same.ss_score)<delta)].shape[0]
    Nbar_SS[i] = data_diff[((SS_scores[i]-data_diff.ss_score)>=0) & ((SS_scores[i]-data_diff.ss_score)<delta)].shape[0]

# determine conditional probability
P_SS_F = N_SS / data_same.shape[0]
P_SS_Fbar = Nbar_SS / data_diff.shape[0]

# Plot values of conditional probabilities
fig, ax = plt.subplots(1,1,figsize=(16,9))
ax.plot(SS_scores, P_SS_F, '^r-', lw=2, alpha = 0.6, label='Same superfamily')
ax.plot(SS_scores, P_SS_Fbar, '2b-', lw=2, alpha = 0.6, label='Different superfamily')
ax.set_xlim([0,1])
ax.legend(loc='best', frameon=False)
plt.show()
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/SS_cond_probs.svg', format="svg")





# # for sampled subset
# N_F = 14140496
# N_Fbar = 373434065

N_F = data_same.shape[0]
N_Fbar = data_diff.shape[0]

P_F = N_F / (N_F + N_Fbar)
P_Fbar = 1 -P_F


SS_scores = np.linspace(0.0,1.0,101)
delta = SS_scores[1]-SS_scores[0]
N_SS = np.zeros_like(SS_scores)
Nbar_SS = np.zeros_like(SS_scores)

# count up all pairs that fall in the interval bins
for i in range(len(SS_scores)):
    N_SS[i] = data_same[((SS_scores[i]-data_same.ss_score)>=0) & ((SS_scores[i]-data_same.ss_score)<delta)].shape[0]
    Nbar_SS[i] = data_diff[((SS_scores[i]-data_diff.ss_score)>=0) & ((SS_scores[i]-data_diff.ss_score)<delta)].shape[0]
# determine conditional probability
P_SS_F = N_SS / data_same.shape[0]
P_SS_Fbar = Nbar_SS / data_diff.shape[0]

P_F_SS = (P_SS_F * P_F) / (P_SS_F * P_F + P_SS_Fbar * P_Fbar)
P_Fbar_SS = (P_SS_Fbar * P_Fbar) / (P_SS_F * P_F + P_SS_Fbar * P_Fbar)


plt.rcParams["axes.labelsize"] = 26
plt.rcParams["xtick.labelsize"] = 20
plt.rcParams["ytick.labelsize"] = 20
plt.rcParams["figure.titlesize"] = 29

# Plot values of posterior probabilities
fig, ax = plt.subplots(1,1,figsize=(16,9))
ax.plot(SS_scores, P_F_SS, '^r-', lw=2, alpha = 0.6, label='Same superfamily')
ax.plot(SS_scores, P_Fbar_SS, 'ob-', lw=2, alpha = 0.6, label='Different superfamily')
ax.set_xlim([0,1])
ax.legend(loc='best', frameon=False)
plt.show()
plt.savefig('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_s20_plots/posteriori_probs_small_numbers.svg', format="svg")




#############


# Create x values from 0 to 100, with interval 20
x = np.linspace(0, 100, 6)

# Create y values using the quadratic function y = x + 0.01 * x^2
y = x + 0.01 * x**2

# Plot the curve
plt.figure(figsize=(8, 8))
plt.plot(x, y, label='Non-linear curve', linewidth=2)

# Add grid, ticks, and labels
plt.grid(True)
plt.xticks(np.arange(0, 101, 20))
plt.yticks(np.arange(0, 101, 20))
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Non-linear Curve from 0 to 100')

# Add legend
plt.legend()

# Show the plot
plt.show()
