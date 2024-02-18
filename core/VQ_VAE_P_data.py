# VQ-VAE 

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from tensorflow import keras
from tensorflow.keras import layers
#import tensorflow_probability as tfp
import tensorflow as tf
# import librosa 
# import librosa.display as ld
import sklearn
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, ParameterGrid
import pickle
from tensorflow.keras.layers import Dropout
from keras.callbacks import EarlyStopping
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
#from statistics import covariance
import multiprocessing as mp
from sklearn.preprocessing import MinMaxScaler
import random



class VectorQuantizer(layers.Layer):
    def __init__(self, num_embeddings, embedding_dim, beta=0.25, **kwargs):
        super().__init__(**kwargs)
        self.embedding_dim = embedding_dim
        self.num_embeddings = num_embeddings
        
        print("num_embeddings,embedding_dim",self.num_embeddings,self.embedding_dim)
        # The `beta` parameter is best kept between [0.25, 2] as per the paper.
        self.beta = beta

        # Initialize the embeddings which we will quantize.
        w_init = tf.random_uniform_initializer()
        self.embeddings = tf.Variable(initial_value=w_init(shape=(self.embedding_dim,self.num_embeddings), dtype="float32"),
            trainable=True,
            name="embeddings_vqvae",
        )
        
    # to save the model     
    def get_config(self):

        config = super().get_config()
        config.update({
            'num_embeddings': self.num_embeddings,
            'embedding_dim': self.embedding_dim,
            'beta': self.beta,
        })
        return config
    
    @classmethod
    def from_config(cls, config):
        return cls(**config)
    
    
    def call(self, x):
        # Calculate the input shape of the inputs and
        # then flatten the inputs keeping `embedding_dim` intact.
        input_shape = tf.shape(x)
        
        print("input_shape........",input_shape)
        
        flattened = tf.reshape(x, [-1, self.embedding_dim])
        
        print("flattened input_shape........",flattened.shape)
        
        
        # Quantization.
        encoding_indices = self.get_code_indices(flattened)
        
        print("encoding_indices.......shape",encoding_indices.shape)
        
        encodings = tf.one_hot(encoding_indices, self.num_embeddings) #minimum distance to the corresponding encoder output 
                                                                      #is mapped as one and the remaining codes are mapped as zeros
        quantized = tf.matmul(encodings, self.embeddings, transpose_b=True)

        # Reshape the quantized values back to the original input shape
        quantized = tf.reshape(quantized, input_shape)


        # Calculate vector quantization loss and add that to the layer. You can learn more
        # about adding losses to different layers here:
       
        commitment_loss = tf.reduce_mean((tf.stop_gradient(quantized) - x) ** 2)
        codebook_loss = tf.reduce_mean((quantized - tf.stop_gradient(x)) ** 2)
        self.add_loss(self.beta * commitment_loss + codebook_loss)

        # Straight-through estimator.
        quantized = x + tf.stop_gradient(quantized - x)
        return quantized

    def get_code_indices(self, flattened_inputs):
        # Calculate L2-normalized distance between the inputs and the codes .
        
        print("flattened_inputs..embeddings shape",flattened_inputs.shape,self.embeddings.shape)
        
        similarity = tf.matmul(flattened_inputs, self.embeddings)
        #similarity = tf.matmul(flattened_inputs, tf.transpose(self.embeddings))
        
        
        #print("similarity.......",similarity,len(similarity))
        
        # for i in self.embeddings:
        #  print("self.embeddings.......",i)

        distances = (
            tf.reduce_sum(flattened_inputs ** 2, axis=1, keepdims=True)
            + tf.reduce_sum(self.embeddings ** 2, axis=0)
            - 2 * similarity
        )
        #print("distances.......",distances)

        # Derive the indices for minimum distances.
        encoding_indices = tf.argmin(distances, axis=1)
        
        return encoding_indices
    


def get_encoder(latent_dim=2):
    
    encoder = tf.keras.models.Sequential(name="encoder")
    encoder_inputs = encoder.add(keras.Input(shape=(4,)))
    encoder.add(layers.Dense(2,activation='relu')) # works better without any activation, although relu is ok.  use_bias=True by default
    encoder.add(Dropout(0.4))  # gives the best results
    encoder_outputs = encoder.add(layers.Dense(2))
    
    
    return encoder





def get_decoder(latent_dim=2):
  
  decoder = tf.keras.models.Sequential(name="decoder")   
  latent_inputs = decoder.add(keras.Input(shape=(2,)))
  decoder.add(layers.Dense(2,activation='relu'))
  decoder.add(Dropout(0.4)) # gives the best results
  decoder_outputs = decoder.add(layers.Dense(4))

  
  return decoder




def get_vqvae(latent_dim=2, num_embeddings=26):
    vq_layer = VectorQuantizer(num_embeddings, latent_dim, name="vector_quantizer")
    encoder = get_encoder(latent_dim)
    decoder = get_decoder(latent_dim)
    inputs = keras.Input(shape=(4,))
    encoder_outputs = encoder(inputs)
    print("get_vqvae .............vq_layer.",vq_layer)
    quantized_latents = vq_layer(encoder_outputs)
    outputs = decoder(quantized_latents)
    
    return keras.Model(inputs, outputs, name="vq_vae")

print(get_vqvae().summary())


class VQVAETrainer(keras.models.Model):
    def __init__(self, train_variance, latent_dim=2, num_embeddings=26, **kwargs):
        super(VQVAETrainer, self).__init__(**kwargs)
        self.train_variance = train_variance
        self.latent_dim = latent_dim
        self.num_embeddings = num_embeddings

        self.vqvae = get_vqvae(self.latent_dim, self.num_embeddings)

        self.total_loss_tracker = keras.metrics.Mean(name="total_loss")
        self.reconstruction_loss_tracker = keras.metrics.Mean(
            name="reconstruction_loss"
        )
        self.vq_loss_tracker = keras.metrics.Mean(name="vq_loss")

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.vq_loss_tracker,
        ]

    def train_step(self, x):
        
        with tf.GradientTape() as tape:
            # Outputs from the VQ-VAE.
            reconstructions = self.vqvae(x)

            # Calculate the losses.
            reconstruction_loss = (
                tf.reduce_mean((x - reconstructions) ** 2) / self.train_variance
            )
            total_loss = reconstruction_loss + sum(self.vqvae.losses)

        # Backpropagation.
        grads = tape.gradient(total_loss, self.vqvae.trainable_variables)
        self.optimizer.apply_gradients(zip(grads, self.vqvae.trainable_variables))

        # Loss tracking.
        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.vq_loss_tracker.update_state(sum(self.vqvae.losses))

        # Log results.
        return {
            "loss": self.total_loss_tracker.result(),
            "reconstruction_loss": self.reconstruction_loss_tracker.result(),
            "vqvae_loss": self.vq_loss_tracker.result(),
        }




def distance(lista, listb):
    return sum( (b - a) ** 2 for a,b in zip(lista, listb) ) ** .5


# if __name__ == "__main__": 

#"cath_final_features_with_BM.pickle"    
with open("cath_final_features_with_BM_without_AA.pickle", "rb") as handle:

        all_proteins_dictionary = pickle.load(handle)

   
###suffle
l = list(all_proteins_dictionary.items())
random.shuffle(l)
all_proteins_dictionary = dict(l)




train_size = int(0.8 * len(all_proteins_dictionary))

train_set_dict = dict(list(all_proteins_dictionary.items())[:train_size])

with open("cath_final_features_with_BM_without_AA_train_set(0.8).pickle", "wb") as handle:
        pickle.dump(train_set_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


test_set_dict = dict(list(all_proteins_dictionary.items())[train_size:])

with open("cath_final_features_with_BM_without_AA_test_set(0.2).pickle", "wb") as handle:
        pickle.dump(test_set_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



print(len(all_proteins_dictionary),len(train_set_dict),len(test_set_dict))


# to return a group of the key-value
# pairs in the dictionary
values = train_set_dict.values()

# Convert object to a list
data = list(values)
    
#inputs = np.array([])
input_list = []
for i in data:
    for j in i:
    
     input_list.append(j)
   

inputs_np = np.array(input_list,dtype="float32")     
print("all inputs_np shape",inputs_np.shape)



# train_size = int(0.8 * len(inputs_np))

# train_set = inputs_np[:train_size]

train_set = inputs_np

# print("train_set 10",train_set[:5])

# train_set_df = pd.DataFrame(train_set, columns = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
#         "THR","TRP","TYR","VAL","H", "S", "L","distance"])

# # print("train_set_df 10",train_set_df.head(10))

# one_hot_column = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
#         "THR","TRP","TYR","VAL","H", "S", "L"]

train_set_df = pd.DataFrame(train_set, columns = ["H", "S", "L","distance"])

# print("train_set_df 10",train_set_df.head(10))

one_hot_column = ["H", "S", "L"]

column_to_scale = ["distance"]

scaler = MinMaxScaler()
train_scaled_columns  = scaler.fit_transform(train_set_df[column_to_scale]) 

# print("scaled_columns 10",scaled_columns[:5])
train_one_hot_encoded_column = train_set_df[one_hot_column]
train_set_s = np.concatenate([train_one_hot_encoded_column,train_scaled_columns], axis=1)

# print("train_set_s 10",train_set_df_s[:5])

#inve_s = scaler.inverse_transform(train_set_s)

# to return a group of the key-value
# pairs in the dictionary
test_values = test_set_dict.values()

# Convert object to a list
test_data = list(test_values)
    
#inputs = np.array([])
test_input_list = []
for i in test_data:
    for j in i:
    
     test_input_list.append(j)
   

test_inputs_np = np.array(test_input_list,dtype="float32")     
print("all test_inputs_np shape",test_inputs_np.shape)

test_set = test_inputs_np

# test_set = inputs_np[train_size:]

# test_set_df = pd.DataFrame(test_set, columns = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
#         "THR","TRP","TYR","VAL","H", "S", "L","distance"])

test_set_df = pd.DataFrame(test_set, columns = ["H", "S", "L","distance"])

test_scaled_columns  = scaler.fit_transform(test_set_df[column_to_scale])
test_one_hot_encoded_column = test_set_df[one_hot_column]

test_set_s = np.concatenate([test_one_hot_encoded_column,test_scaled_columns], axis=1)

print(".....train_set.......len(test_set)",len(test_set))


data_variance = np.var(train_set)
#data_variance = np.var(train_set_s)


#Train the VQ-VAE model
es = EarlyStopping(monitor='loss', mode='min', verbose=1, patience=5)
vqvae_trainer = VQVAETrainer(data_variance, latent_dim=2, num_embeddings=26)
vqvae_trainer.compile(optimizer=keras.optimizers.Adam(learning_rate=0.001)) # Defaults 0.001 works best
# history = vqvae_trainer.fit(train_set, epochs=500, batch_size=32,callbacks=[es]) # batch_size 32 gives best result
history = vqvae_trainer.fit(train_set_s, epochs=30, batch_size=64,callbacks=[es]) # batch_size 32 gives best result


# # Reconstruction results on the test set 
trained_vqvae_model = vqvae_trainer.vqvae
print("trained_vqvae_model",type(trained_vqvae_model))

#save the model
trained_vqvae_model.save('vqvae_model_train_.8_percent_4-2_c-2-26_scaled')      

