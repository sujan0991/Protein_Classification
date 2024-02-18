import os
import subprocess
import itertools




file_seq = '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/cath_non_redundant_pdbs_s40_sequence_ForkPoolWorker-1.csv'
path_tostore='/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_all_TM_after_ck1-2-3/'
path_pdbs= '/group/bioinf_protstr/Ballal/dompdb_cath_s40/'





PATH_PATTERNS = [
    # working directory
    lambda pdb_id: os.path.join(".", f"{pdb_id}.cif"),
    lambda pdb_id: os.path.join(".", f"{pdb_id}.pdb"),
    # general pdb dir (e.g. for mounting via singularity)
    lambda pdb_id: os.path.join("/pdb", pdb_id),
    lambda pdb_id: os.path.join("/pdb", f"{pdb_id}.pdb"),
    lambda pdb_id: os.path.join("/pdb", f"{pdb_id}.cif"),
    lambda path: path
]



def get_file_path(pdb_id):
    for pattern in PATH_PATTERNS:
        path = pattern(pdb_id)
        if os.path.exists(path):
            return path
        
US_WIN = os.path.join( "/group/bioinf_protstr/Ballal/multiple-protein-structure-alignment/mpsa/us_align/", "USAlign.exe")
US = os.path.join( "/group/bioinf_protstr/Ballal/multiple-protein-structure-alignment/mpsa/us_align/", "USAlign")


def get_row_path(row):
    return get_file_path(row["pdb_id"])

def get_cigar(seq_a, seq_b, aln):
    cigar = ""
    for index in range(len(aln)):
        a = seq_a[index:index + 1]
        b = seq_b[index:index + 1]
        q = aln[index:index + 1]
        if a == "-":
            cigar += "D"
            continue
        if b == "-":
            cigar += "I"
            continue
        if q == ".":
            cigar += "M"
        if q == ":":
            cigar += "P"
    return "".join([str(len(list(group))) + l for l, group in itertools.groupby(cigar)])


def run_aln_cmd(command, row_source, row_target, fast=False):
    try:
        fp_source = row_source
        fp_target = row_target
    except FileNotFoundError as e:
        return int(0), float(1e6), float(0), float(0), str("")
    command.append(fp_source)
    command.append(fp_target)
    command.append("-outfmt")
    command.append("0")
    if fast:
        command.append("-fast")
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()
    p.terminate()
    if err:
        return int(0), float(1e6), float(0), float(0), str("")
    output = output.decode("ascii")
    output = output.replace(r"\r", "").split("\n")
    l_1 = output[10].split(" ")[-2]
    l_2 = output[11].split(" ")[-2]
    aliLen, rmsd, seq_id = output[13][:-1].split(",")
    aliLen, rmsd, seq_id = aliLen.split("=")[-1].strip(), rmsd.split("=")[-1].strip(), seq_id.split("=")[-1].strip()
    tm_1 = output[14].split(" ")[1]
    tm_2 = output[15].split(" ")[1]
    seq_a, aln, seq_b = output[-6:-3]
    cigar = get_cigar(seq_a, seq_b, aln)
    #print(int(aliLen), float(rmsd), float(tm_1), float(tm_2), str(cigar))
    return int(aliLen), float(rmsd), float(tm_1), float(tm_2), str(cigar)




def align_us_(row_source, row_target):
    command = [US_WIN if os.name == "nt" else US]
    return run_aln_cmd(command, row_source, row_target)



#align_us_("/group/bioinf_protstr/Ballal/cath_trimmed_PDBs2/4zutA01.cif", "/group/bioinf_protstr/Ballal/cath_trimmed_PDBs2/1uatA00.cif")










############


import pandas as pd
import sys
import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np
import multiprocessing
import sys
from Levenshtein import distance as lev







#####
b= int(sys.argv[1]) ##### 0......9
no_threads = int(sys.argv[2])
no_batches=5

import warnings
warnings.filterwarnings("ignore")



#file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_ForkPoolWorker-2.csv'

df_seq= pd.read_csv(file_seq,skipinitialspace=True)
df_seq =df_seq[['id']]
P2= set(df_seq.index)


def main(group_indices):
    current_process=multiprocessing.current_process().name
    counter = 0
    for i in group_indices:
        counter = counter + 1 
        df_p1= df_seq.loc[i] # first protein
        p_name=df_p1.id
        ### load pdb p_name
        P2_sub_index = [x for x in  P2 if x >i]
        df_p2= df_seq[df_seq.index.isin(P2_sub_index)] # keep only the proteins of interest to compare with
        df_p2['TM_query']=-9
        df_p2['TM_subject']=-9
        df_p2['query_P']=p_name
        print(current_process,counter,len(group_indices), len(df_p2))
        sys.stdout.flush()
        for j in df_p2.index:
            p2_name=df_p2.loc[j,'id']
            p1=path_pdbs+p_name+'.pdb'
            p2=path_pdbs+p2_name+'.pdb'
            try:
                align=align_us_(p1,p2)
                df_p2.loc[j,'TM_query'] =  align[2]
                df_p2.loc[j,'TM_subject'] =   align[3]
                # df_p2.loc[j,'TM_min'] =  min(align[2],align[3])
                # df_p2.loc[j,'TM_max'] =   max(align[2],align[3])
            except:
                df_p2.loc[j,'TM_query'] =-8
                df_p2.loc[j,'TM_subject'] =-8
            #print(p_name,p2_name, align[2:4])
        df_p2 = df_p2[['id', 'query_P','TM_query','TM_subject']]
        df_p2.to_csv(path_tostore+'queryP_'+p_name +'batch_'+ str(b)+'.csv', index=False)
        


pool = multiprocessing.Pool(no_threads)   
P1= df_seq.index.to_list()
groups_indicesGlobal = np.array_split(P1, no_batches)  # int(no_batches/no_threads)
groups_indices = np.array_split(groups_indicesGlobal[b], no_threads)  # int(no_batches/no_threads)
pool.map(main, groups_indices)


exit ()

## cat *csv > /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_s40_results.csv


import pandas as pd
df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_s40_results.csv')
df = df.query("`query_P` != 'query_P'")
df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/CATH_TM_s40_results.csv')

