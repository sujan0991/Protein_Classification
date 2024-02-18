import numpy as np
import pandas as pd
import glob
import os


def compress_string(s):
    compressed = s[0]
    for char in s[1:]:
        if char != compressed[-1]:
            compressed += char
    return compressed



## 
df_seq= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/ForkPoolWorker-1.csv')
df_seq=df_seq[~df_seq.SS_seq.isna()]
# df_24k = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/domain_id_list_24k.csv')
# list_24k = df_24k['domain_id'].tolist()
# df_seq = df_seq[df_seq['id'].isin(list_24k)]
# df2 = df_seq.query('AA_seq.str.contains("\?")')
#df2_ids = df2['id'].tolist()
#df_seq = df_seq[~df_seq['id'].isin(df2_ids)]
df_seq['SS_seq_reduced'] = -9
for index, row in df_seq.iterrows():
    input_string = row['SS_seq']
    compressed_result = compress_string(input_string)
    df_seq.loc[index,'SS_seq_reduced'] = compressed_result
    print(df_seq)

#df_seq.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_seq_reduced_ForkPoolWorker-1.csv',index=False)
df_seq.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/Rad52_motif_sequence_full_reduced.csv',index=False)



####  for AF
exit()

chunk_folder = "/group/bioinf_protstr/All_alphaFoldSS_Jun2023/codes/chunks_alphaFold/data_chunks/"

ck_id=14000 ## must change this 


file_list = os.listdir('/group/bioinf_protstr/All_alphaFoldSS_Jun2023/codes/chunks_alphaFold/data_chunks/')
   
df = pd.DataFrame(file_list, columns=["id"])

for index, row in df.iloc[14000:].iterrows():
    csv_name = row['id']
    chunk_file = chunk_folder + csv_name
    chunk_df = pd.read_csv(chunk_file)
    chunk_df.columns=['id','SS_seq']
    chunk_df['SS_seq_reduced'] = -9
    ck_id+=1
    #print(chunk_df)
    for index, row in chunk_df.iterrows():
        input_string = row['SS_seq']
        compressed_result = compress_string(input_string)
        chunk_df.loc[index,'SS_seq_reduced'] = compressed_result
    chunk_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/AF_sequences_full_reduced_v2/chunk_'+ str(ck_id)+'.csv',index=False)




