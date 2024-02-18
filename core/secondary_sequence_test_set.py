import pandas as pd
import sys

import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np



def convert_to_string(s):
 
    # initialization of string to ""
    new = ""
 
    # traverse in the string
    for x in s:
        new += x
 
    # return string
    return new


##"cath_all_proteins_dictionary_with_BM_speed_test.pickle"
## "domain_sequence_SS_speed_test.csv"

with open("cath_all_proteins_dictionary_with_BM_speed_test.pickle", "rb") as handle:

        all_proteins_dictionary = pickle.load(handle)




# lev_df2 = pd.read_csv('distance_comparison_test_set(0.2)_10k_superfamily_domain_without_AA.csv')


####i don't know why it's different
# test_set_min_4_S_F = pd.read_csv("test_set_2_10k_superfamily_domain.csv")
    
# print(list(test_set_min_4_S_F.columns.values))


# keysList = lev_df2.domain1.unique().tolist()


# print("key list len",len(keysList))


# original_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
# key_list = ['a', 'c']
# new_dict = {k: original_dict[k] for k in key_list if k in original_dict}
# print(new_dict)
# # Output: {'a': 1, 'c': 3}

# test_set_dict = {k: all_proteins_dictionary[k] for k in keysList if k in all_proteins_dictionary}

# ss_seq_dict = dict()

# for k in keysList:

#         # GET THE values from all_proteins_dictionary
#         val = all_proteins_dictionary[k]
        
#         #print("ss seq",val["ss_struc"])
        
#         ss_seq = convert_to_string(val["ss_struc"])
        
#         ss_seq_dict[k] = ss_seq
        

# print("ss seq dict len",len(ss_seq_dict))        
            
# df = pd.DataFrame(list(ss_seq_dict.items()), columns=['domain_id','sequence'])   

# df.to_csv('sequence_original_ss_test_set_10k_superfamily_domain.csv',index=False)      


####### all sequence

ss_seq_dict = dict()       

for k, v in all_proteins_dictionary.items():

    ss_seq = convert_to_string(v["ss_struc"])
    ss_seq_dict[k] = ss_seq     
    
print("ss seq dictt len",len(ss_seq_dict))    

            
df = pd.DataFrame(list(ss_seq_dict.items()), columns=['domain_id','sequence'])   

df.to_csv('domain_sequence_SS_speed_test.csv',index=False)      

