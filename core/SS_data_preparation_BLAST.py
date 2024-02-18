import pandas as pd
import sys
import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np
import multiprocessing
import sys
import glob
import itertools
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq




df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_ss_seq_superfamily_after_ck1-2-3.csv')

## S > X
## L>Y
## H>Z

df['SS_seq'] = df['SS_seq'].str.replace('S','X')
df['SS_seq'] = df['SS_seq'].str.replace('L','Y')
df['SS_seq'] = df['SS_seq'].str.replace('H','Z')


## create fasta

seq_recs = []

for id,seq in df[['id','SS_seq']].to_records(index=False):
    if  seq.isalpha():
        seq_ = Seq(seq)
        seq_ =  SeqIO.SeqRecord(seq_, id, description = "")
        seq_recs.append(seq_)
    else:
        print('contain other charecter...........................',id,seq)



    
SeqIO.write(seq_recs,'/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_ss_seq_superfamily_after_replace_S_H_L.fasta', 'fasta')



## seperate fasta 

## cd /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/all_single_SS_fasta_file
## awk '/^>/{filename=substr($0,2) ".fasta"}; {print > filename}' cath_ss_seq_superfamily_after_replace_S_H_L.fasta