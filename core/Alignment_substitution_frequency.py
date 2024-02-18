import pandas as pd
import matplotlib.pyplot as plt
import pickle 

def count_substitution_frequency(string1, string2, results):
    if len(string1) != len(string2):
        print("The strings must be of equal length for alignment.")
        return
    for i in range(len(string1)):
        if string1[i] == 'H' and string2[i] == 'S':
            results['H_to_S'] += 1
        elif string1[i] == 'S' and string2[i] == 'H':
            results['H_to_S'] += 1  
        elif string1[i] == 'H' and string2[i] == 'L':
            results['H_to_L'] += 1 
        elif string1[i] == 'L' and string2[i] == 'H':
            results['H_to_L'] += 1 
        elif string1[i] == 'S' and string2[i] == 'L':
            results['S_to_L'] += 1 
        elif string1[i] == 'L' and string2[i] == 'S':
            results['S_to_L'] += 1 
        elif string1[i] == 'H' and string2[i] == '-':
            results['H_to_B'] += 1 
        elif string1[i] == '-' and string2[i] == 'H':
            results['H_to_B'] += 1
        elif string1[i] == 'S' and string2[i] == '-':
            results['S_to_B'] += 1 
        elif string1[i] == '-' and string2[i] == 'S':
            results['S_to_B'] += 1 
        elif string1[i] == 'L' and string2[i] == '-':
            results['L_to_B'] += 1 
        elif string1[i] == '-' and string2[i] == 'L':
            results['L_to_B'] += 1                        


df = pd.read_csv('/group/bioinf/Ballal/cath_s20_local_results_full_reduced.csv',
                 usecols=['local_ident_full','query_aln_string_full','target_aln_string_full'])

results = {'H_to_S': 0, 'H_to_L': 0, 'S_to_L': 0, 
           'H_to_B': 0,  'S_to_B': 0,  'L_to_B': 0}

for row in df.itertuples(index=False):
    string1 = row.query_aln_string_full
    string2 = row.target_aln_string_full
    count_substitution_frequency(string1, string2, results)

with open('/group/bioinf/Ballal/substitution_frequency_local_custom_s20_all.pkl', 'wb') as f:
    pickle.dump(results, f)

with open('/group/bioinf/Ballal/substitution_frequency_local_custom_scop_all.pkl', 'rb') as f:
    loaded_dict = pickle.load(f)

results = {'H_to_S': 0, 'H_to_L': 0, 'S_to_L': 0, 
           'H_to_B': 0,  'S_to_B': 0,  'L_to_B': 0}

results['H_to_S'] = loaded_dict['H_to_S'] + loaded_dict['S_to_H']
results['H_to_L'] = loaded_dict['H_to_L'] + loaded_dict['L_to_H']
results['S_to_L'] = loaded_dict['S_to_L'] + loaded_dict['L_to_S']
results['H_to_B'] = loaded_dict['H_to_B'] + loaded_dict['B_to_H']
results['S_to_B'] = loaded_dict['S_to_B'] + loaded_dict['B_to_S']
results['L_to_B'] = loaded_dict['L_to_B'] + loaded_dict['B_to_L']

names = list(results.keys())
values = list(results.values())

plt.figure(figsize=(20,6))
plt.bar(range(len(results)), values, tick_label=names)
plt.savefig('/group/bioinf/Ballal/plots/sub_frequency_local_custom_scop_all.svg', format="svg")
plt.close()

exit()

# Example 
string1 = "-------------------------------------------------------------------------------------------------------------LLLLL--HHHHHHHHHHHHHHHHHHHHHHHHHHHHLLLLLLLLSSSSLLLLLSSSSSLLLLLLLLLLLLLLSSSSSSLLSSSSSL"
string2 = "LLLLHHHHHHHLLLLLLLLLLLLLLHHHHHLLLLLHHHHHHHHHHHHLLLLLHHHHHHHHHHHHHHHHLHHHHHHHHHHHHHLLLHHHHHHHHHHHHLHHHHHLLHHHHLLLLLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHLLHHHLLLHHHHHLLHHHHHHHHHHHHHHHHHHHHHHHHHLLLL------"
count_substitution_frequency(string1, string2, results)


S1 = 'LHHHHHHHHHH----------HHHHHHHHHHHHHLSSSLLLLLSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSLLLLLHHHHHHHHHHLLLLSSSSSSSSSSSLLLLLLLLLL-LLLSSSSSSSSSSLLLLLLLHHHHHHHHHHHHHHHHHHHHHHHHHLLLLLLLLLLLSSSSHHHHHHHL--LLLLHHHHHHHHHHHHLSSSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSLLLLSSSSSSSSSSLLLHHHHHHHHHHHLLHHHHHLHHHHHHHHLLLLLSSSSSSSHHHHHHHHHLLLLHHHHLLLLLLHHHHHHLLLLL'
S2 = 'L----------LSSSSSSLLLHHHHHHHHHHHHH------LLLSSS--------------------------LSSS---LLHHHHHHHHHH-------------------------HLLL-----SSSSSLLLLLLL---------------HHHHHHHHHH------LLLLLSSSS--------SSLLLL--HHHHHHHHHH------------------------------------------------LLL------SSSS-LL-----------LLHHHHH-HHHHHHHH-LLLL----------------LLLL--HH------HHHHHH----L'
count_substitution_frequency(S1, S2,results)





