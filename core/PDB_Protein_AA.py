import pandas as pd
import sys

import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np




with open("cath_all_proteins_dictionary_with_BM.pickle", "rb") as handle:

        all_proteins_dictionary = pickle.load(handle)

   
 ######### convert to one hot vector #########

    # main 20 amino acids
keys_aminos = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ] 


final_features = dict()

for k in all_proteins_dictionary.keys():

        # GET THE values from all_proteins_dictionary
        val = all_proteins_dictionary[k]
        # to store features in one hot format
        features = []

        for k2 in range(len(val["resns"])):

            featuresDic_aminos = dict.fromkeys(
                keys_aminos, 0
            )  # dictionary of amino acids name, the form of one hot vector.First, set all 20 bits(20 amino acid) to 0
            
            # if the resns matches the name in keys_aminos, set that bit to 1, rest will be 0
            if val["resns"][k2] in keys_aminos:

                featuresDic_aminos[val["resns"][k2]] = 1

               

                # append all three features in features
                features.append(
                    list(featuresDic_aminos.values())
                    
                )

            # if not,should we store them in 21st bit?????
            else:

                print("have a different amino acid", k, val["resns"][k2])
                

                # we r not taking this aa in features

        # store all feature in final_features where the key is the PDB id and value id features list
        final_features[k] = features

    
with open("cath_final_features_AA_with_BM.pickle", "wb") as handle:
        pickle.dump(final_features, handle, protocol=pickle.HIGHEST_PROTOCOL)
