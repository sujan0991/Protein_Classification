import pandas as pd
import numpy as np
import pickle

# Open the pickle file and load the data
with open("cath_all_proteins_dictionary_SS_AA_with_BM.pickle", "rb") as handle:
    all_proteins_dictionary = pickle.load(handle)

# Initialize an empty dictionary for amino acid sequences
aa_seq_dict = {}

# Iterate over the items in the dictionary
for k, v in all_proteins_dictionary.items():
    seq = v['resns']
    print('seq...', seq)
    
    # Process the sequence
    seq = seq.strip('\n')
    seq = seq.replace("\n", "")
    seq = seq.split("_")
    seq = seq[2]
    seq = seq[1:]
    print('seq...????', seq)
    
    aa_seq_dict[k] = seq

print(len(aa_seq_dict))

# Create a DataFrame from the dictionary
df = pd.DataFrame(list(aa_seq_dict.items()), columns=['domain_id', 'sequence'])

print(df.head(), len(df))

# Write the DataFrame to a CSV file
df.to_csv("domain_sequence_AA_safe_cath.csv", index_label="index")
