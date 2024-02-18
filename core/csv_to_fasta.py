from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


#df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv')
df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_scope_ForkPoolWorker-1.csv')

df.drop_duplicates(subset='id', keep="last", inplace=True)
df=df[~df.SS_seq.isna()]
df=df[~df.AA_seq.isna()]
df=df[~df.id.isna()]
df.reset_index(inplace=True)

print("....",len(df))

# df_24k = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/domain_id_list_24k.csv')
# list_24k = df_24k['domain_id'].tolist()

# df = df[df['id'].isin(list_24k)]



# df2 = df.query('AA_seq.str.contains("\?")')
# df2_ids = df2['id'].tolist()

# df = df[~df['id'].isin(df2_ids)]

df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/cath_seq_random_1000_speed_test.csv')

seq_recs = []



for id,seq in df[['id','AA_seq']].to_records(index=False):
    if  seq.isalpha():
        seq_ = Seq(seq)
        seq_ =  SeqIO.SeqRecord(seq_, id, description = "")
        seq_recs.append(seq_)
    else:
        print('contain other charecter...........................',id,seq)

print(len(seq_recs))##s40:31846     , scop: 9138


    
SeqIO.write(seq_recs,'/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/speed_test_for_1000_aa_seq.fasta', 'fasta')


    
    
   