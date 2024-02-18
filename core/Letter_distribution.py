import pandas as pd
import pickle
import time
import numpy as np
import itertools 
import matplotlib.pyplot as plt
import seaborn as sns
import math
from collections import Counter


with open("cath_all_proteins_dictionary_with_BM.pickle", "rb") as handle:

        all_proteins_dictionary = pickle.load(handle)
        


df = pd.read_csv("domain_sequence_without_AA_10k_test_set(0.2).csv")      


print(len(all_proteins_dictionary),len(df))  

lt = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

helix_A = 0
sheet_A = 0
loop_A = 0

helix_B = 0
sheet_B = 0
loop_B = 0

helix_C = 0
sheet_C = 0
loop_C = 0

helix_D = 0
sheet_D = 0
loop_D = 0

helix_E = 0
sheet_E = 0
loop_E = 0

helix_F = 0
sheet_F = 0
loop_F = 0

helix_G = 0
sheet_G = 0
loop_G = 0

helix_H = 0
sheet_H = 0
loop_H = 0

helix_I = 0
sheet_I = 0
loop_I = 0

helix_J = 0
sheet_J = 0
loop_J = 0

helix_K = 0
sheet_K = 0
loop_K = 0

helix_L = 0
sheet_L = 0
loop_L = 0

helix_M = 0
sheet_M = 0
loop_M = 0


helix_N = 0
sheet_N = 0
loop_N = 0

helix_O = 0
sheet_O = 0
loop_O = 0

helix_P = 0
sheet_P = 0
loop_P = 0

helix_Q = 0
sheet_Q = 0
loop_Q = 0

helix_R = 0
sheet_R = 0
loop_R = 0

helix_S = 0
sheet_S = 0
loop_S = 0

helix_T = 0
sheet_T = 0
loop_T = 0

helix_U = 0
sheet_U = 0
loop_U = 0

helix_V = 0
sheet_V = 0
loop_V = 0

helix_W = 0
sheet_W = 0
loop_W = 0

helix_X = 0
sheet_X = 0
loop_X = 0

helix_Y = 0
sheet_Y = 0
loop_Y = 0

helix_Z = 0
sheet_Z = 0
loop_Z = 0

distance_var_dict = {}
#aa_var_dict = {}

for single_lt in lt:
    
    distance_var_dict[f"{single_lt}"] = []
    #aa_var_dict[f"aa_{single_lt}"] = []
    

for k, v in all_proteins_dictionary.items():

       #print("result..........i",k, type(v))
       
        seq = df.loc[df['domain_id'] == k,'sequence']
        
        seq_string = np.array_str(seq.values)
        seq_string = seq_string[2:]
        seq_string = seq_string[:-2]
        
        
        # aa = v['resns']
        ss = v['ss_struc']
        distance = v['distances']
        
        print(len(seq_string))
        
        seq_unique = set(seq_string)
        
        print(seq_unique)
        
        
        for i,v in enumerate(distance):
            
            if len(seq_string) == len(distance):
                
            #   print(seq_string[i],v,type(aa[i]))
              
              if seq_string[i] == 'A':
                  
                  distance_var_dict['A'] += [v]
                  #aa_var_dict["aa_A"] += [aa[i]]
                 
              elif seq_string[i] == 'B':
                  
                  distance_var_dict['B'] += [v]
                  #aa_var_dict["aa_B"] += [aa[i]]
               
              elif seq_string[i] == 'C':
                  
                  distance_var_dict['C'] += [v]
                  #aa_var_dict["aa_C"] += [aa[i]]
                  
              elif seq_string[i] == 'D':
                  
                  distance_var_dict['D'] += [v]
                  #aa_var_dict["aa_D"] += [aa[i]]
                  
              elif seq_string[i] == 'E':
                  
                  distance_var_dict['E'] += [v]
                  #aa_var_dict["aa_E"] += [aa[i]]
                  
              elif seq_string[i] == 'F':
                  
                  distance_var_dict['F'] += [v]
                  #aa_var_dict["aa_F"] += [aa[i]]
                  
              elif seq_string[i] == 'G':
                  
                  distance_var_dict['G'] += [v]
                  #aa_var_dict["aa_G"] += [aa[i]]
                  
              elif seq_string[i] == 'H':
                  
                  distance_var_dict['H'] += [v]
                  #aa_var_dict["aa_H"] += [aa[i]]
                  
              elif seq_string[i] == 'I':
                  
                  distance_var_dict['I'] += [v]
                  #aa_var_dict["aa_I"] += [aa[i]]
                  
              elif seq_string[i] == 'J':
                  
                  distance_var_dict['J'] += [v] 
                  #aa_var_dict["aa_J"] += [aa[i]] 
                                                
              elif seq_string[i] == 'K':
                  
                  distance_var_dict['K'] += [v]
                  #aa_var_dict["aa_K"] += [aa[i]]
                  
              elif seq_string[i] == 'L':
                  
                  distance_var_dict['L'] += [v]
                  #aa_var_dict["aa_L"] += [aa[i]]
                  
              elif seq_string[i] == 'M':
                  
                  distance_var_dict['M'] += [v]
                  #aa_var_dict["aa_M"] += [aa[i]]
                  
              elif seq_string[i] == 'N':
                  
                  distance_var_dict['N'] += [v]
                  #aa_var_dict["aa_N"] += [aa[i]]
                  
              elif seq_string[i] == 'O':
                  
                  distance_var_dict['O'] += [v]
                  #aa_var_dict["aa_O"] += [aa[i]]
                  
              elif seq_string[i] == 'P':
                  
                  distance_var_dict['P'] += [v]
                  #aa_var_dict["aa_P"] += [aa[i]]
                  
              elif seq_string[i] == 'Q':
                  
                  distance_var_dict['Q'] += [v]
                  #aa_var_dict["aa_Q"] += [aa[i]]
                  
              elif seq_string[i] == 'R':
                  
                  distance_var_dict['R'] += [v]
                  #aa_var_dict["aa_R"] += [aa[i]]
                  
              elif seq_string[i] == 'S':
                  
                  distance_var_dict['S'] += [v]
                  #aa_var_dict["aa_S"] += [aa[i]]
                  
              elif seq_string[i] == 'T':
                  
                  distance_var_dict['T'] += [v]
                  #aa_var_dict["aa_T"] += [aa[i]]
                  
              elif seq_string[i] == 'U':
                  
                  distance_var_dict['U'] += [v]
                  #aa_var_dict["aa_U"] += [aa[i]]
                  
              elif seq_string[i] == 'V':
                  
                  distance_var_dict['V'] += [v]
                  #aa_var_dict["aa_V"] += [aa[i]]
                  
              elif seq_string[i] == 'W':
                  
                  distance_var_dict['W'] += [v]
                  #aa_var_dict["aa_W"] += [aa[i]]
                  
              elif seq_string[i] == 'X':
                  
                  distance_var_dict['X'] += [v]
                  #aa_var_dict["aa_X"] += [aa[i]]
                  
              elif seq_string[i] == 'Y':
                  
                  distance_var_dict['Y'] += [v]
                  #aa_var_dict["aa_Y"] += [aa[i]]
                  
              elif seq_string[i] == 'Z':
                  
                  distance_var_dict['Z'] += [v]  
                  #aa_var_dict["aa_Z"] += [aa[i]]                          
                                                   
#print("..........",aa_var_dict)     

  

        


# exit()

# ############## ss ###############
        
        for i,v in enumerate(ss):
            
            if len(seq_string) == len(ss):
                
              
              
              if seq_string[i] == 'A' and v == 'H' :
                  
                 helix_A = helix_A + 1 
                 
              elif seq_string[i] == 'A' and v == 'S' :   
                  
                  sheet_A = sheet_A + 1
                  
              elif seq_string[i] == 'A' and v == 'L' :   
                  
                  loop_A = loop_A + 1    
                  
                  
                  
              elif seq_string[i] == 'B' and v == 'H' :   
                  
                  helix_B = helix_B + 1
                  
              elif seq_string[i] == 'B' and v == 'S' :   
                  
                  sheet_B = sheet_B + 1    
              elif seq_string[i] == 'B' and v == 'L' :   
                  
                  loop_B = loop_B + 1
              
              
                  
              elif seq_string[i] == 'C' and v == 'H' :   
                  
                  helix_C = helix_C + 1
                  
              elif seq_string[i] == 'C' and v == 'S' :   
                  
                  sheet_C = sheet_C + 1    
              elif seq_string[i] == 'C' and v == 'L' :   
                  
                  loop_C = loop_C + 1    
                  
                  
              
              elif seq_string[i] == 'D' and v == 'H' :   
                  
                  helix_D = helix_D + 1
                  
              elif seq_string[i] == 'D' and v == 'S' :   
                  
                  sheet_D = sheet_D + 1    
              elif seq_string[i] == 'D' and v == 'L' :   
                  
                  loop_D = loop_D + 1
                  
                  
              
              elif seq_string[i] == 'E' and v == 'H' :   
                  
                  helix_E = helix_E + 1
                  
              elif seq_string[i] == 'E' and v == 'S' :   
                  
                  sheet_E = sheet_E + 1    
              elif seq_string[i] == 'E' and v == 'L' :   
                  
                  loop_E = loop_E + 1
                  
                  
              elif seq_string[i] == 'F' and v == 'H' :   
                  
                  helix_F = helix_F + 1
                  
              elif seq_string[i] == 'F' and v == 'S' :   
                  
                  sheet_F = sheet_F + 1    
              elif seq_string[i] == 'F' and v == 'L' :   
                  
                  loop_F = loop_F + 1
              
              
              elif seq_string[i] == 'G' and v == 'H' :   
                  
                  helix_G = helix_G + 1
                  
              elif seq_string[i] == 'G' and v == 'S' :   
                  
                  sheet_G = sheet_G + 1    
              elif seq_string[i] == 'G' and v == 'L' :   
                  
                  loop_G = loop_G + 1
                  
                  
              elif seq_string[i] == 'H' and v == 'H' :   
                  
                  helix_H = helix_H + 1
                  
              elif seq_string[i] == 'H' and v == 'S' :   
                  
                  sheet_H = sheet_H + 1    
              elif seq_string[i] == 'H' and v == 'L' :   
                  
                  loop_H = loop_H + 1
                  
                  
              elif seq_string[i] == 'I' and v == 'H' :   
                  
                  helix_I = helix_I + 1
                  
              elif seq_string[i] == 'I' and v == 'S' :   
                  
                  sheet_I = sheet_I + 1    
              elif seq_string[i] == 'I' and v == 'L' :   
                  
                  loop_I = loop_I + 1
                  
                  
                  
              elif seq_string[i] == 'J' and v == 'H' :   
                  
                  helix_J = helix_J + 1
                  
              elif seq_string[i] == 'J' and v == 'S' :   
                  
                  sheet_J = sheet_J + 1    
              elif seq_string[i] == 'J' and v == 'L' :   
                  
                  loop_J = loop_J + 1                        
             
             
              elif seq_string[i] == 'K' and v == 'H' :   
                  
                  helix_K = helix_K + 1
                  print("helix_K..............")
              elif seq_string[i] == 'K' and v == 'S' :   
                  print("sheet_K..............????")
                  sheet_K = sheet_K + 1    
              elif seq_string[i] == 'K' and v == 'L' :   
                  
                  loop_K = loop_K + 1  
                  
                  
              elif seq_string[i] == 'L' and v == 'H' :   
                  
                  helix_L = helix_L + 1
                  
              elif seq_string[i] == 'L' and v == 'S' :   
                  
                  sheet_L = sheet_L + 1    
              elif seq_string[i] == 'L' and v == 'L' :   
                  
                  loop_L = loop_L + 1  
                  
                  
                  
              elif seq_string[i] == 'M' and v == 'H' :   
                  
                  helix_M = helix_M + 1
                  
              elif seq_string[i] == 'M' and v == 'S' :   
                  
                  sheet_M = sheet_M + 1    
              elif seq_string[i] == 'M' and v == 'L' :   
                  
                  loop_M = loop_M + 1  
                  
                  
                  
              elif seq_string[i] == 'N' and v == 'H' :   
                  
                  helix_N = helix_N + 1
                  
              elif seq_string[i] == 'N' and v == 'S' :   
                  
                  sheet_N = sheet_N + 1    
                  
              elif seq_string[i] == 'N' and v == 'L' :   
                  
                  loop_N = loop_N + 1  
                  
                  
                  
              elif seq_string[i] == 'O' and v == 'H' :   
                  
                  helix_O = helix_O + 1
                  
              elif seq_string[i] == 'O' and v == 'S' :   
                  
                  sheet_O = sheet_O + 1    
              elif seq_string[i] == 'O' and v == 'L' :   
                  
                  loop_O = loop_O + 1  
                  
                  
                  
              elif seq_string[i] == 'P' and v == 'H' :   
                  
                  helix_P = helix_P + 1
                  
              elif seq_string[i] == 'P' and v == 'S' :   
                  
                  sheet_P = sheet_P + 1    
              elif seq_string[i] == 'P' and v == 'L' :   
                  
                  loop_P = loop_P + 1  
                  
                  
              elif seq_string[i] == 'Q' and v == 'H' :   
                  
                  helix_Q = helix_Q + 1
                  
              elif seq_string[i] == 'Q' and v == 'S' :   
                  
                  sheet_Q = sheet_Q + 1    
              elif seq_string[i] == 'Q' and v == 'L' :   
                  
                  loop_Q = loop_Q + 1  
                  
                  
              elif seq_string[i] == 'R' and v == 'H' :   
                  
                  helix_R = helix_R + 1
                  
              elif seq_string[i] == 'R' and v == 'S' :   
                  
                  sheet_R = sheet_R + 1    
              elif seq_string[i] == 'R' and v == 'L' :   
                  
                  loop_R = loop_R + 1  
                  
                  
              elif seq_string[i] == 'S' and v == 'H' :   
                  
                  helix_S = helix_S + 1
                  
              elif seq_string[i] == 'S' and v == 'S' :   
                  
                  sheet_S = sheet_S + 1    
              elif seq_string[i] == 'S' and v == 'L' :   
                  
                  loop_S = loop_S + 1     
                  
                  
              elif seq_string[i] == 'T' and v == 'H' :   
                  
                  helix_T = helix_T + 1
                  
              elif seq_string[i] == 'T' and v == 'S' :   
                  
                  sheet_T = sheet_T + 1    
              elif seq_string[i] == 'T' and v == 'L' :   
                  
                  loop_T = loop_T + 1   
                  
                  
                  
              elif seq_string[i] == 'U' and v == 'H' :   
                  
                  helix_U = helix_U + 1
                  
              elif seq_string[i] == 'U' and v == 'S' :   
                  
                  sheet_U = sheet_U + 1    
              elif seq_string[i] == 'U' and v == 'L' :   
                  
                  loop_U = loop_U + 1   
                  
                  
              elif seq_string[i] == 'V' and v == 'H' :   
                  
                  helix_V = helix_V + 1
                  
              elif seq_string[i] == 'V' and v == 'S' :   
                  
                  sheet_V = sheet_V + 1    
              elif seq_string[i] == 'V' and v == 'L' :   
                  
                  loop_V = loop_V + 1   
                  
                  
                  
              elif seq_string[i] == 'W' and v == 'H' :   
                  
                  helix_W = helix_W + 1
                  
              elif seq_string[i] == 'W' and v == 'S' :   
                  
                  sheet_W = sheet_W + 1    
              elif seq_string[i] == 'W' and v == 'L' :   
                  
                  loop_W = loop_W + 1   
                  
                  
              elif seq_string[i] == 'X' and v == 'H' :   
                  
                  helix_X = helix_X + 1
                  
              elif seq_string[i] == 'X' and v == 'S' :   
                  
                  sheet_X = sheet_X + 1    
              elif seq_string[i] == 'X' and v == 'L' :   
                  
                  loop_X = loop_X + 1   
                  
                  
                  
              elif seq_string[i] == 'Y' and v == 'H' :   
                  
                  helix_Y = helix_Y + 1
                  
              elif seq_string[i] == 'Y' and v == 'S' :   
                  
                  sheet_Y = sheet_Y + 1    
              elif seq_string[i] == 'Y' and v == 'L' :   
                  
                  loop_Y = loop_Y + 1   
                  
                  
                  
              elif seq_string[i] == 'Z' and v == 'H' :   
                  
                  helix_Z = helix_Z + 1
                  
              elif seq_string[i] == 'Z' and v == 'S' :   
                  
                  sheet_Z = sheet_Z + 1    
              elif seq_string[i] == 'Z' and v == 'L' :   
                  
                  loop_Z = loop_Z + 1                                                            
             
            else:
             print('......................',i)     
        



helix_dict = {'A':helix_A,'B':helix_B,'C':helix_C,'D':helix_D,'E':helix_E,'F':helix_F,'G':helix_G,'H':helix_H,'I':helix_I,'J':helix_J,
              'K':helix_K,'L':helix_L,'M':helix_M,'N':helix_N,'O':helix_O,'P':helix_P,'Q':helix_Q,'R':helix_R,'S':helix_S,'T':helix_T,
              'U':helix_U,'V':helix_V,'W':helix_W,'X':helix_X,'Y':helix_Y,'Z':helix_Z}

sheet_dict = {'A':sheet_A,'B':sheet_B,'C':sheet_C,'D':sheet_D,'E':sheet_E,'F':sheet_F,'G':sheet_G,'H':sheet_H,'I':sheet_I,'J':sheet_J,
              'K':sheet_K,'L':sheet_L,'M':sheet_M,'N':sheet_N,'O':sheet_O,'P':sheet_P,'Q':sheet_Q,'R':sheet_R,'S':sheet_S,'T':sheet_T,
              'U':sheet_U,'V':sheet_V,'W':sheet_W,'X':sheet_X,'Y':sheet_Y,'Z':sheet_Z}

loop_dict = {'A':loop_A,'B':loop_B,'C':loop_C,'D':loop_D,'E':loop_E,'F':loop_F,'G':loop_G,'H':loop_H,'I':loop_I,'J':loop_J,
              'K':loop_K,'L':loop_L,'M':loop_M,'N':loop_N,'O':loop_O,'P':loop_P,'Q':loop_Q,'R':loop_R,'S':loop_S,'T':loop_T,
              'U':loop_U,'V':loop_V,'W':loop_W,'X':loop_X,'Y':loop_Y,'Z':loop_Z}



# print('helix_A.........',helix_dict,sheet_dict,loop_dict)    



with open("helix_letter.pickle", "wb") as handle:
        pickle.dump(helix_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
with open("sheet_letter.pickle", "wb") as handle:
        pickle.dump(sheet_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
with open("loop_letter.pickle", "wb") as handle:
        pickle.dump(loop_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        



###### plot ########    
       
with open("distance_letter.pickle", "wb") as handle:
        pickle.dump(distance_var_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
# # with open("aa_letter.pickle", "wb") as handle:
# #         pickle.dump(aa_var_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)              
  
        
        
with open("distance_letter.pickle", "rb") as handle:

        distance_letter = pickle.load(handle) 
        

#print("distance_letter....",distance_letter.items())        

# for key,n_key in zip(distance_letter.keys(), lt):
#     distance_letter[n_key] = distance_letter.pop(key)
    

labels, data = [*zip(*distance_letter.items())]  # 'transpose' items to parallel key, value lists

# or backwards compatable    
labels, data = distance_letter.keys(), distance_letter.values()

plt.boxplot(data)
plt.xticks(range(1, len(labels) + 1), labels)
plt.savefig('distance_letter.png')
       
       
       
       
       
       
        
with open("helix_letter.pickle", "rb") as handle:

        helix = pickle.load(handle)        
        
aa=([np.log10(i) for i in helix.values()])         

fig, ax = plt.subplots(2)


ax[0].bar(list(helix.keys()), aa)
    #plt.savefig('aa_count.png')


ax[1].bar(list(helix.keys()), helix.values())
plt.savefig('Helix_letter.png')
#plt.show()        



with open("sheet_letter.pickle", "rb") as handle:

        helix = pickle.load(handle)        
        
aa=([np.log10(i) for i in helix.values()])         

fig, ax = plt.subplots(2)


ax[0].bar(list(helix.keys()), aa)
    #plt.savefig('aa_count.png')


ax[1].bar(list(helix.keys()), helix.values())
plt.savefig('sheet_letter.png')




with open("loop_letter.pickle", "rb") as handle:

        helix = pickle.load(handle)        
        
aa=([np.log10(i) for i in helix.values()])         

fig, ax = plt.subplots(2)


ax[0].bar(list(helix.keys()), aa)
    #plt.savefig('aa_count.png')


ax[1].bar(list(helix.keys()), helix.values())
plt.savefig('loop_letter.png')



