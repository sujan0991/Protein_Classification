import glob
import itertools
import glob
import pandas as pd
import os
import numpy as np



#### check2 

id_list = []
outcome_of_ck2_list = []

path= '/group/bioinf_protstr/ali/cifCheck2/' # use your path
all_files = glob.glob((path+'*.txt'))
for filename in all_files:
     id = filename.split("cifCheck2/")[1][:4]
     id_list.append(id)
     ## check if 1st row is '_pdbx_poly_seq_scheme.hetero '
     with open(filename) as f:
      first_line = f.readline().strip('\n')
      if first_line == '_pdbx_poly_seq_scheme.hetero ':  
          df = pd.DataFrame([line.strip().split() for line in f.readlines()])
          column_5_value = df.iloc[0,5]
          column_6_value = df.iloc[0,6]
          if column_5_value == column_6_value and column_6_value != '?':
              #print('no_conflict')
              outcome_of_ck2_list.append('no_conflict')
          else:
              #print('conflict')
              outcome_of_ck2_list.append('conflict')    
      else:
          print('noData')  
          outcome_of_ck2_list.append('noData')
          
result_df = pd.DataFrame(list(zip(id_list,outcome_of_ck2_list)), columns=['PDB_Id','outcome_of_ck2'])          
result_df.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/pdb_check_2_list.csv')

###########


####### check3
id_list = []
outcome_of_ck3_list = []

no_conflict = 0
conflict = 0
noData = 0

path= '/group/bioinf_protstr/ali/cifCheck2/' # use your path
all_files = glob.glob((path+'*.txt'))
for filename in all_files:
     id = filename.split("cifCheck2/")[1][:4]
     id_list.append(id)
     ## check if 1st row is '_pdbx_poly_seq_scheme.hetero '
     with open(filename) as f:
      first_line = f.readline().strip('\n')
      if first_line == '_pdbx_poly_seq_scheme.hetero ':  
          df = pd.DataFrame([line.strip().split() for line in f.readlines()])
          if df.iloc[len(df)-1,0] == '#':
            df = df[:-1]
            df['ck3'] = np.where(((df[5]==df[6]) & (df[5] != '?') & (df[6] != '?')), 'confirmed', 'not_confirmed')
            if df['ck3'].eq('confirmed').all():
              no_conflict += 1
              outcome_of_ck3_list.append('no_conflict')
            else:
              outcome_of_ck3_list.append('conflict')    
              conflict += 1
          else:
            df['ck3'] = np.where(((df[5]==df[6]) & (df[5] != '?') & (df[6] != '?')), 'confirmed', 'not_confirmed')
            if df['ck3'].eq('confirmed').all():
              no_conflict += 1
              outcome_of_ck3_list.append('no_conflict')
            else:
              outcome_of_ck3_list.append('conflict')    
              conflict += 1 
      else:
          print('noData')  
          outcome_of_ck3_list.append('noData')
          noData += 1
          
 
result_df = pd.DataFrame(list(zip(id_list,outcome_of_ck3_list)), columns=['PDB_Id','ck3'])          
result_df.to_csv('/group/bioinf_protstr/ali/ck3.csv')
 





######## merge ck3 and ck2-ck1

ck1_ck2_df =  pd.read_csv('/group/bioinf_protstr/ali/ck1-ck2_merged.csv')
ck3_df =  pd.read_csv('/group/bioinf_protstr/ali/ck3.csv')
 
ck1_ck2_df= pd.merge(ck1_ck2_df, ck3_df ,how='left', on="PDB_Id") 
ck1_ck2_df = ck1_ck2_df[['PDB_Id', 'ck2', 'ck1', 'ck3']]
ck1_ck2_df.groupby(['ck2','ck1','ck3']).count()
                    
ck1_ck2_df.to_csv('/group/bioinf_protstr/ali/ck1-ck2-ck3_merged.csv')
 
           
           

                    
####### safe pdb

ck_df = pd.read_csv('/group/bioinf_protstr/ali/ck1-ck2-ck3_merged.csv')

idx = np.where((ck_df['ck3']=='no_conflict') & (ck_df['ck2']=='no_conflict') & (ck_df['ck1']=='no_conflict'))
no_con_no_con_df = ck_df.loc[idx]
no_con_no_con_df.to_csv('/group/bioinf_protstr/ali/ck1-ck2-ck3_no_conflict.csv')

idx = np.where((ck_df['ck3']=='no_conflict') & (ck_df['ck2']=='no_conflict') & (ck_df['ck1']=='noData'))
no_con_no_data_df = ck_df.loc[idx]
no_con_no_data_df.to_csv('/group/bioinf_protstr/ali/ck2-ck3_no_conflict_ck1_noData.csv')
