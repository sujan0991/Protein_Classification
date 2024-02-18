from time import process_time_ns
import numpy as np
import tensorflow as tf
from tensorflow import keras
import pandas as pd 
import pickle
#import VQ_VAE_Test_P_data_2 as vq

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import pearsonr

from Levenshtein import distance as lev
from itertools import combinations
import multiprocessing as mp
import itertools
from sklearn.preprocessing import MinMaxScaler




#"cath_final_features_with_BM.pickle"   
with open("cath_final_features_with_BM.pickle", "rb") as handle:

        all_proteins_dictionary = pickle.load(handle)
    

    #vqvae_model_18-12_c-12-26_scaled
model = tf.keras.models.load_model("vqvae_model_18-12_c-12-26_scaled",compile=False)    


print(model.summary())

quantizer = model.get_layer("vector_quantizer")

codebook = quantizer.embeddings

print("codebook......??.",codebook.shape)

codebook2 = tf.transpose(codebook)
print("codebook2......??.",codebook2.shape)



lt = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
lt_dict = {}
for i, s_VQ in enumerate(codebook2):
    
    vq = tf.get_static_value(s_VQ)
    lt_dict[lt[i]] = vq
    



encoder = model.get_layer("encoder")

#'11baA00', '11baB00', '11bgA00', '11bgB00'  11gsA01



def get_encoded_indecies(domain_id):
    
    print("get_encoded_indecies")
    
    domain_value = list(all_proteins_dictionary[domain_id])
    
    scaler = MinMaxScaler()
     # fit and transform in one step
    domain_value_s = scaler.fit_transform(domain_value)
     
    
    encoded_outputs = encoder.predict(domain_value_s)
    flat_enc_outputs = encoded_outputs.reshape(-1, encoded_outputs.shape[-1])
    #print("flat_enc_outputs shape",flat_enc_outputs.shape)
    similarity = tf.matmul(flat_enc_outputs, codebook)
    
    distances = (
            tf.reduce_sum(flat_enc_outputs ** 2, axis=1, keepdims=True)
            + tf.reduce_sum(codebook ** 2, axis=0)
            - 2 * similarity)
    
    encoding_indices = tf.argmin(distances, axis=1)
    
    return  encoding_indices
   
def get_distance(indies1,indies2):
    
    print("get_sequence")
    
    sequence = ''


    for i in indies1:
   
      for j, s_VQ in enumerate(codebook2):
         if i == j:
          
            vq = tf.get_static_value(s_VQ)
            
            for k, v in lt_dict.items():
                
                if (v == vq).all():
                    sequence = sequence + k
                    
    sequence2 = ''


    for i in indies2:
   
      for j, s_VQ in enumerate(codebook2):
         if i == j:
          
            vq = tf.get_static_value(s_VQ)
            
            for k, v in lt_dict.items():
                
                if (v == vq).all():
                    sequence2 = sequence2 + k           
                    
    distance = lev(sequence, sequence2)   
    
    return distance                                      


###############


test_p_dict =  list(all_proteins_dictionary.keys())[:500]

result_list = []

print("result_list",type(result_list))

all_comb = itertools.combinations(test_p_dict,2)

for domain_key in all_comb:
    
    
    print("domain_key......",domain_key)
    domain1_key = domain_key[0]
    domain2_key = domain_key[1]
    indies1 = get_encoded_indecies(domain1_key)
    indies2 = get_encoded_indecies(domain2_key)
    
    dist = get_distance(indies1,indies2)
    
    print("distance between",domain1_key,domain2_key,"is", dist)
    
    return_tuple = (domain1_key,domain2_key,dist)
    
    result_list.append(return_tuple)

df = pd.DataFrame(result_list, columns=['id','dimain2', 'distance'])
 
df2 = pd.read_csv("cath-domain-description-v4_3_0.csv")
  
merge_result = pd.merge(df, df2, on="id")
 
 #change the order of column
merge_result = merge_result[['id','Class','Arch','Topol','Homol','distance','dimain2']]
 
 #rename column
merge_result.columns = ['domain1','Class','Arch','Topol','Homol','distance','id']
 
#  merge_result = merge_result.drop(['index'], axis=1)  
merge_result.to_csv("distance_comparison.csv", index_label="index")
 
df3 = pd.read_csv("distance_comparison.csv")
  
full_result = pd.merge(df3, df2, on="id")             
        

full_result = full_result.drop(['index_x'], axis=1)
full_result = full_result.drop(['index_y'], axis=1) 
  #rename column
full_result.columns = ['domain1','Class1','Arch1','Topol1','Homol1','lev_distance','domain2','Class2','Arch2','Topol2','Homol2']
 
full_result.to_csv("lev_distance_comparison_500.csv", index_label="index")

###############

# def get_domain_key(domain_key):
    
#     print("domain_key......",domain_key)
#     domain1_key = domain_key[0]
#     domain2_key = domain_key[1]
#     indies1 = get_encoded_indecies(domain1_key)
#     indies2 = get_encoded_indecies(domain2_key)
    
#     dist = get_distance(indies1,indies2)
    
#     print("distance between",domain1_key,domain2_key,"is", dist)
#     return domain1_key,domain2_key, dist


# sliced_p_dict =  list(all_proteins_dictionary.keys())[:10000]

# print("sliced_p_dict.......",len(sliced_p_dict))
# #test_p_dict =  ['11baA00', '11baB00', '11bgA00', '11bgB00','11gsA01', '11gsA02', '11gsB01', '11gsB02']


# if __name__ == "__main__":
    
#  print("inside main")   
      
#  with mp.Pool(10) as pool:

#         # result = pool.map(get_domain_key, itertools.combinations(all_proteins_dictionary,2)) 
#         result = pool.map(get_domain_key, itertools.combinations(sliced_p_dict,2)) 
        
        
#         #print("result......",result)
        
        
#  df = pd.DataFrame(result, columns=['id','dimain2', 'distance'])
 
#  df2 = pd.read_csv("cath-domain-description-v4_3_0.csv")
  
#  merge_result = pd.merge(df, df2, on="id")
 
#  #change the order of column
#  merge_result = merge_result[['id','Class','Arch','Topol','Homol','distance','dimain2']]
 
#  #rename column
#  merge_result.columns = ['domain1','Class','Arch','Topol','Homol','distance','id']
 
# #  merge_result = merge_result.drop(['index'], axis=1)  
#  merge_result.to_csv("distance_comparison.csv", index_label="index")
 
#  df3 = pd.read_csv("distance_comparison.csv")
  
#  full_result = pd.merge(df3, df2, on="id")             
        

#  full_result = full_result.drop(['index_x'], axis=1)
#  full_result = full_result.drop(['index_y'], axis=1) 
#   #rename column
#  full_result.columns = ['domain1','Class1','Arch1','Topol1','Homol1','lev_distance','domain2','Class2','Arch2','Topol2','Homol2']
 
#  full_result.to_csv("lev_distance_comparison_50k.csv", index_label="index")
 
 
 
    

    
    
 