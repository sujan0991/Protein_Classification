#from time import process_time_ns
import numpy as np
# import tensorflow as tf
# from tensorflow import keras
import pandas as pd 
from Levenshtein import distance as lev
from itertools import combinations
import multiprocessing as mp
import itertools
import time



def get_distance(domain_key):
    domain1_key = domain_key[0]
    domain2_key = domain_key[1]
    print(domain1_key,domain2_key)
    sequence1 = df.loc[df.id == domain1_key,'SS_seq'].item()
    sequence2 = df.loc[df.id == domain2_key,'SS_seq'].item()
    print(sequence1,',',sequence2)
    ss_distance = lev(sequence1, sequence2)
    #distance = levenshtein(sequence1, sequence2)   
    seq1_len = len(sequence1)
    seq2_len = len(sequence2)
    print("......",ss_distance)
    return domain1_key,domain2_key, seq1_len, seq2_len,ss_distance





# read and concate dataframe
# df = pd.read_csv("sequence_test_10k_superfamily_domain.csv") # actual number 9028

#df = pd.read_csv("domain_sequence_SS_speed_test.csv") 
file_seq= '/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_child_ForkPoolWorker-2.csv'
df= pd.read_csv(file_seq)
df=df[~df.SS_seq.isna()]
df=df[~df.AA_seq.isna()]

df.drop_duplicates(subset='id', keep="last", inplace=True)


df =df.sample(1000)
#domain_sequence_without_AA_10k_test_set.csv

print("domain_sequence_SS_speed_test.csv.csv len",len(df))



id = df["id"]

test_p_list =  list(id)

print("test_p_list length",len(test_p_list))




if __name__ == "__main__":
    
 start_time = time.time()
   
 print("inside main")
      
 
    
 with mp.Pool(1) as pool:
        result = pool.map(get_distance, itertools.combinations(test_p_list,2)) 
        
        
 print("result......",len(result))
 print("lev distance calculation time.....",time.time() - start_time)
 
 #domain1_key,domain2_key, seq1_len, seq2_len,ss_distance,aa_distance
 distance_df = pd.DataFrame(result, columns=['id','domain2','sequence1_len','sequence2_len','ss_distance'])
 #df = pd.DataFrame(result, columns=['id','sequence1','dimain2','sequence2','distance'])
 distance_df.to_csv("distance_comparison_SS_safe_cath.csv", index_label="index")
 
 
