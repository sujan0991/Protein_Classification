
#@markdown Installing PyMOL takes a while. Only needed for image generation.
install_pymol = True #@param {type:"boolean"}
install_biotite = True #@param {type:"boolean"}
install_biopython = True #@param {type:"boolean"}


##############################################################################################
#
# Imports

from multiprocessing import Pool, Manager
import pandas as pd
import os,time
from statistics import median,mean,stdev    # mean and median AlphaFold confidence
import matplotlib.pyplot as plt             # plot histograms etc.
import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pymol.cmd as cmd
import numpy as np
import glob
from multiprocessing import Pool
from Levenshtein import distance as lev
import psutil

##############################################################################################

# Set the number of threads
num_threads = 30
scores = (2,-1,-3,-2)

# Paths to query and target PDB file folders
folder_path_target = "/group/bioinf/Users/Stella/lab_project/af/"
folder_path_query = "/group/bioinf_protstr/Ballal/Research_Project/Stella_project/motif_pdb/"





# Function to process a single PDB file and return the AA,SS string
def process_pdb_file(pdb_path):
    try:
        pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
        cmd.load(pdb_path, pdb_id)
        seq = cmd.get_fastastr(pdb_id)
        sec_struc = []
        selection = "n. CA"
        cmd.iterate(selection, "sec_struc.append(ss)", space=locals())
        ss_string = ''.join(['L' if a == '' else a for a in sec_struc])
        return pdb_id, ss_string,seq
    except Exception as e:
        print(f"Error processing PDB file: {str(e)}", file=sys.stderr)
        return None
    finally:
        cmd.delete(pdb_id)
        sys.stdout.flush()

# Function to perform sequence alignment and return the alignment score
def perform_alignment(args):
    Q_string, T_string = args
    alignment_score = lev(Q_string, T_string, weights=(1,2,1))
    longer_seq_length = max(len(Q_string), len(T_string))
    long_score = (1 - alignment_score / longer_seq_length)

    shorter_seq_length = min(len(Q_string), len(T_string))
    subtract_len = abs(len(Q_string) - len(T_string))
    short_score = (1 - (alignment_score- subtract_len)/ shorter_seq_length)

    score_fit = long_score

    if len(Q_string) > len(T_string):
        score_Q = long_score
        score_T = short_score
    else:
        score_Q = short_score
        score_T = long_score
    return score_Q,score_T,score_fit

# def perform_alignment(Q_string, T_string ,c):
#     alignment_score_levG = lev(Q_string, T_string, weights=(abs(c[0]), abs(c[1]), abs(c[2])))
#     levG_ident = alignment_score_levG / max(len(Q_string), len(T_string))
#     local_alignments = pairwise2.align.localms(Q_string, T_string, c[0], c[1], c[2], c[3])
#     if len(local_alignments) > 0:
#         best_local_alignment = max(local_alignments, key=lambda alignment: alignment[2])
#     else:
#         best_local_alignment = [0, 0, 0]
#     best_local_alignment_ident = best_local_alignment[2] / max(len(Q_string), len(T_string))

#     global_alignments = pairwise2.align.globalms(Q_string, T_string, c[0], c[1], c[2], c[3])
#     best_global_alignment = max(global_alignments, key=lambda alignment: alignment[2])

#     return best_local_alignment[2], best_local_alignment_ident, levG_ident



# Shared list to store the alignment scores
manager = Manager()
alignment_scores = manager.list()

# Function to process a single PDB file and perform alignment with all target PDB files
def process_file(f_T):
    if f_T.endswith(".pdb"):
        pdb_path_Target = os.path.join(folder_path_target, f_T)
        pdb_id_T, ss_string_T, aa_string_T = process_pdb_file(pdb_path_Target)
        ss_string_T = ss_string_T.replace('B', 'L').replace('N', 'L')
        #print(f"File: {f_T}")
        #print(f"PDB ID: {pdb_id_T}")
        #print(f"Secondary Structure: {ss_string_T}")
        #print()

        # Iterate over each QUERY PDB file
        for f_Q in os.listdir(folder_path_query):
            if f_Q.endswith(".pdb"):
                pdb_path_Query = os.path.join(folder_path_query, f_Q)
                pdb_id_Q, ss_string_Q , aa_string_Q= process_pdb_file(pdb_path_Query)
                ss_string_Q = ss_string_Q.replace('B', 'L').replace('N', 'L')
                #print(f"File: {f_Q}")
                #print(f"PDB ID: {pdb_id_Q}")
                #print(f"Secondary Structure: {ss_string_Q}")

                # Perform sequence alignment and store the alignment score
                alignment_SS_score_Q,alignment_SS_score_T,alignment_SS_score_fit = perform_alignment((ss_string_Q,ss_string_T))
                alignment_AA_score_Q,alignment_AA_score_T,alignment_AA_score_fit = perform_alignment((aa_string_Q,aa_string_T))

  
                #print(f"Alignment Score: {alignment_score}")
                #print()
                alignment_scores.append([str(pdb_id_Q)+'-ssQ', pdb_id_T, alignment_SS_score_Q])
                alignment_scores.append([str(pdb_id_Q)+'-ssT', pdb_id_T,alignment_SS_score_T])
                alignment_scores.append([str(pdb_id_Q)+'-ssF', pdb_id_T,alignment_SS_score_fit])

                alignment_scores.append([str(pdb_id_Q)+'-aaQ', pdb_id_T, alignment_AA_score_Q])
                alignment_scores.append([str(pdb_id_Q)+'-aaT', pdb_id_T,alignment_AA_score_T])
                alignment_scores.append([str(pdb_id_Q)+'-aaF', pdb_id_T,alignment_AA_score_fit])

                

                


# Get the list of target PDB files
target_files = [f for f in os.listdir(folder_path_target) if f.endswith(".pdb")]

# Create a pool of worker processes
pool = Pool(processes=num_threads)

start_time = time.time()
start_cpu_time = time.process_time()

# Process the PDB files in parallel
pool.map(process_file, target_files)


end_time = time.time()
end_cpu_time = time.process_time()
elapsed_time = end_time - start_time
cpu_time = end_cpu_time - start_cpu_time

# Calculate the average speed of alignment
num_sequences = len(alignment_scores)
average_speed = elapsed_time / num_sequences

# Calculate the CPU usage
cpu_usage = psutil.cpu_percent()

print("\nAverage Speed of Alignment: {:.2f} seconds per sequence".format(average_speed))
print("CPU Time: {:.2f} seconds".format(cpu_time))
print("CPU Usage: {:.2f}%".format(cpu_usage))

# Convert the shared list to a regular list
alignment_scores = list(alignment_scores)
alignment_scores = pd.DataFrame(alignment_scores)
alignment_scores.columns = ['Query_Id','Target_Id','Score']
#print(alignment_scores)

n_values = alignment_scores[alignment_scores['Score'] < 0.0]

print(n_values,len(n_values))

#alignment_scores.to_csv('/group/bioinf_protstr/toAli/alignment_scores_radRed.csv')

# Pivot the table to create new columns from unique values in the "Query_Id" column
pivoted_df = alignment_scores.pivot(index="Target_Id", columns="Query_Id", values="Score").reset_index()
#pivoted_df.to_csv('/group/bioinf_protstr/toAli/pivoted_alignment_scores_radRed.csv')



df_info=pd.read_csv('/group/bioinf_protstr/To_ballal/RadRed_allInfo.csv', usecols=['Accession','Superkingdom','Family_name'])
pivoted_df=pivoted_df.merge(df_info, left_on='Target_Id', right_on='Accession', how='inner')
pivoted_df.to_csv('/group/bioinf_protstr/toAli/pivoted_df_alignment_scores_radRed_annotated.csv')



exit()



## best combinations for pairwise in (2,-1,-3,-2), in function it's the return value of best_local_alignment_ident/alignment_SS_score_T
## best lev (ss) is default (1,1,1)
# 









import plotly as plt
import numpy as np
import plotly.express as px
import scipy.stats



df_info=pd.read_csv('/group/bioinf_protstr/To_ballal/RadRed_allInfo.csv', usecols=['Accession','Superkingdom','Family_name'])
pivoted_df=pivoted_df.merge(df_info, left_on='Target_Id', right_on='Accession', how='inner')
pivoted_df.to_csv('/group/bioinf_protstr/Ballal/ss_tm_penalty_test/pivoted_dfalignment_scores_radRed_annotated2.csv-'+str([scores])+'.csv')

df_tm= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/Stella_project/tm_query_subject_5_protein.csv')

df_all = pivoted_df.merge(df_tm[['Unnamed: 0','TMQ_rad52_full','TMS_rad52_full']], left_on='Accession', right_on='Unnamed: 0', how='inner')
df_all['mintm']=df_all[['TMQ_rad52_full','TMS_rad52_full']].min(axis=1)
df_all['maxtm']=df_all[['TMQ_rad52_full','TMS_rad52_full']].max(axis=1)


##corelation betwen ssf,sst,ssq vs tm min, it should be (-)

print('w1,w2,w3',scores)

print('spearmanr and pearsonr tm min vs ssF',scipy.stats.spearmanr(df_all['mintm'],df_all['rad52_motif-ssF']),scipy.stats.pearsonr(df_all['mintm'],df_all['rad52_motif-ssF']))
print('spearmanr and pearsonr tm min vs ssQ',scipy.stats.spearmanr(df_all['mintm'],df_all['rad52_motif-ssQ']),scipy.stats.pearsonr(df_all['mintm'],df_all['rad52_motif-ssQ']))
print('spearmanr and pearsonr tm min vs ssT',scipy.stats.spearmanr(df_all['mintm'],df_all['rad52_motif-ssT']),scipy.stats.pearsonr(df_all['mintm'],df_all['rad52_motif-ssT']))

print('spearmanr and pearsonr tm min vs ssF',scipy.stats.spearmanr(df_all['maxtm'],df_all['rad52_motif-ssF']),scipy.stats.pearsonr(df_all['maxtm'],df_all['rad52_motif-ssF']))
print('spearmanr and pearsonr tm min vs ssQ',scipy.stats.spearmanr(df_all['maxtm'],df_all['rad52_motif-ssQ']),scipy.stats.pearsonr(df_all['maxtm'],df_all['rad52_motif-ssQ']))
print('spearmanr and pearsonr tm min vs ssT',scipy.stats.spearmanr(df_all['maxtm'],df_all['rad52_motif-ssT']),scipy.stats.pearsonr(df_all['maxtm'],df_all['rad52_motif-ssT']))

df_all['color'] = df_all['Family_name']+df_all['Superkingdom']
color=np.array(df_all['color'])
hover_name = df_all.Accession

fig = px.scatter(df_all, x= 'rad52_motif-ssQ',y = 'mintm',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pQ-min_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssQ',y = 'maxtm',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pQ-max_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssQ',y = 'TMQ_rad52_full',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pQ-Q_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssQ',y = 'TMS_rad52_full',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pQ-S_'+str(scores)+'.html', auto_open=False)

##ssq vs tmmax 


fig = px.scatter(df_all, x= 'rad52_motif-ssF',y = 'mintm',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pF-min_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssF',y = 'maxtm',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pF-max_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssF',y = 'TMQ_rad52_full',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pF-Q_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssF',y = 'TMS_rad52_full',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pF-S_'+str(scores)+'.html', auto_open=False)



fig = px.scatter(df_all, x= 'rad52_motif-ssT',y = 'mintm',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pT-min_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssT',y = 'maxtm',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pT-max_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssT',y = 'TMQ_rad52_full',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pT-Q_'+str(scores)+'.html', auto_open=False)
fig = px.scatter(df_all, x= 'rad52_motif-ssT',y = 'TMS_rad52_full',color= color,hover_name=hover_name)
plt.offline.plot(fig, filename = '/group/bioinf_protstr/Ballal/ss_tm_penalty_test/new_plots/pT-S_'+str(scores)+'.html', auto_open=False)






# # Count the number of members from each family overall
# family_counts_all = pivoted_df['Family_name'].value_counts().reset_index()
# family_counts_all.columns = ['Family_name', 'Family_Member_Count_overall']

# # Filter the aligned sequences DataFrame based on the Percentage Score condition
# filtered_df1 = pivoted_df[pivoted_df['Score'] > 0.5]

# # Count the number of members from each family that satisfy the condition
# family_counts_condition = filtered_df1['Family_name'].value_counts().reset_index()
# family_counts_condition.columns = ['Family_name', 'Family_Member_Count_condition']

# # Merge the two DataFrames on the 'family_text' column
# merged_counts = pd.merge(family_counts_all, family_counts_condition, on='Family_name', how='left')
# merged_counts['ratio']= merged_counts.Family_Member_Count_condition/ merged_counts.Family_Member_Count_overall
# merged_counts=merged_counts.replace(np.NaN,0)
# merged_counts.sort_values('ratio',ascending=False, inplace=True)
# # Print the merged DataFrame
# print("Number of members from each family:")
# print(merged_counts)



