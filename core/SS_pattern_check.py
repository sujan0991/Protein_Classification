import pandas as pd





df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/fetched_sequences/df_merged_all_safe_cath_seq_after_ck1-2-3_ForkPoolWorker-1.csv')

df2 = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/domain_id_list_24k.csv')
list_24k = df2['domain_id'].tolist()
df = df[df['id'].isin(list_24k)]

##[]

df['First'] = df['SS_seq'].astype(str).str[0]
df['Last'] = df['SS_seq'].astype(str).str[-1]

df = df.query("`First` != 'L'")
df = df.query("`Last` != 'L'")

df = df[['id','SS_seq']]

lev_tm_df = pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath_lev_tm_HH_blast_blast_ss_results.csv')
lev_tm_df = lev_tm_df[['domain1', 'domain2', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max']]
lev_tm_df.columns = ['id', 'domain2', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max']
lev_tm_df = pd.merge(lev_tm_df, df ,how='left', on="id")
lev_tm_df=lev_tm_df[~lev_tm_df.SS_seq.isna()]
lev_tm_df.columns = ['domain1', 'id', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max', 'SS_seq1']
lev_tm_df = pd.merge(lev_tm_df, df ,how='left', on="id")
lev_tm_df=lev_tm_df[~lev_tm_df.SS_seq.isna()]
lev_tm_df.columns = ['domain1', 'domain2', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max', 'SS_seq1','SS_seq2']
lev_tm_df = lev_tm_df[['domain1', 'domain2', 'cath_superFamily', 'SS_score', 'key_id', 'TM_min', 'TM_max']]





same = lev_tm_df.query('cath_superFamily == 1 & TM_max >= 0.5')
print('>=0.5',len(same))

same = lev_tm_df.query('cath_superFamily == 1 & TM_max < 0.5')
print('< 0.5',len(same))



same = lev_tm_df.query('cath_superFamily == 0 & TM_max >= 0.5')
print('>=0.5',len(same))

same = lev_tm_df.query('cath_superFamily == 0 & TM_max < 0.5')
print('< 0.5',len(same))


