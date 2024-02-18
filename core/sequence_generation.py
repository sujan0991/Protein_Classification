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
import random



#"cath_final_features_with_BM.pickle"   
with open("cath_final_features_with_BM_without_AA.pickle", "rb") as handle:

        all_proteins_dictionary = pickle.load(handle)
    


print("all_proteins_dictionary len",len(all_proteins_dictionary))





    #vqvae_model_18-12_c-12-26_scaled
model = tf.keras.models.load_model("vqvae_model_train_.8_percent_4-2_c-2-26_scaled",compile=False)    







print(model.summary())

quantizer = model.get_layer("vector_quantizer")

codebook = quantizer.embeddings

print("codebook......??.",codebook.shape)

codebook2 = tf.transpose(codebook)
print("codebook2......??.",codebook2.shape)


lt = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
# lt = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e',
# 'e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
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




def get_Sequence(indies):
    
    print("get_sequence")
    
    sequence = ''
    
    for i in indies:
   
      for j, s_VQ in enumerate(codebook2):
         if i == j:
          
            vq = tf.get_static_value(s_VQ)
            
            for k, v in lt_dict.items():
                
                if (v == vq).all():
                    sequence = sequence + k
                    
                    
                    
                    
    #print("se.........",sequence)      
    
    return sequence           
                    
                    
                    
                    
                    
domains_10k = pd.read_csv("test_set_10k_superfamily_domain(0.2).csv")

print(len(domains_10k))
print("",domains_10k.columns.tolist())  
                      
                    
# test_p_dict =  list(all_proteins_dictionary.keys())
test_p_dict = domains_10k['id'].tolist()

# test_p_dict = random.sample(list(all_proteins_dictionary.keys()), 1000)

print("test_p_dict len..........",len(test_p_dict))

result_dict = {}

        

for domain_key in test_p_dict:
    
    
    
    indies = get_encoded_indecies(domain_key)     
    
    seq = get_Sequence(indies)   
    
    
    print("sequence",seq)
    
    result_dict[domain_key] = seq
    
print("result_dict len..........",len(result_dict))
    
df = pd.DataFrame(list(result_dict.items()),columns = ['domain_id','sequence'])

df.to_csv("domain_sequence_without_AA_10k_test_set(0.2).csv", index_label="index")


