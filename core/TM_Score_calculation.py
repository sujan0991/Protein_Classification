import warnings
import tmscoring
warnings.filterwarnings('ignore')
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


b=19


def main(p_combinations_child):
    current_process=multiprocessing.current_process().name
    p_error_list=[]
    counter = 0
    path= '/group/bioinf_protstr/Ballal/cath_trimmed_PDBs/'
    for i in p_combinations_child.index:
        counter = counter + 1 
        if (counter % 1000) == 0:print(current_process,counter,len(p_combinations_child))
        sys.stdout.flush()
        try:
             p1=path+p_combinations_child.loc[i,'id1']+'.pdb'
             p2=path+p_combinations_child.loc[i,'id2']+'.pdb'
             #p1_name=p_combinations_child.loc[i,'id1']
             #p2_name=p_combinations_child.loc[i,'id2']
             alignment1 = tmscoring.TMscoring(p1,p2)
             alignment1.optimise()
             s1= alignment1.tmscore(**alignment1.get_current_values())
             alignment2 = tmscoring.TMscoring(p2,p1)
             alignment2.optimise()
             s2= alignment2.tmscore(**alignment2.get_current_values())
             p_combinations_child.loc[i,'TM_min']=(min(s1,s2))
             p_combinations_child.loc[i,'TM_max']=(max(s1,s2))
        except: 
             p_error_list.append([p_combinations_child.loc[i,'id1'],p_combinations_child.loc[i,'id2']])
    print(current_process, len(p_combinations_child), len(p_error_list))
    pd.DataFrame(p_error_list).to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_childrenFilter_False/'+'/error_list_'+current_process+'.txt', index=False)
    p_combinations_child.to_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations_childrenFilter_False/'+str(b)+'/'+current_process+'.csv', index=False)


print('batch: ', b)
no_threads = 30
no_batches=20
pool = multiprocessing.Pool(no_threads)
p_combinations= pd.read_csv('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/p_combinations35K_secondFilter:False.csv', usecols=['id1','id2'])
p_combinations['TM_min']=-9.1
p_combinations['TM_max']=-9.1
pool = multiprocessing.Pool(no_threads)

groupsGlobal = np.array_split(p_combinations, no_batches)  # int(no_batches/no_threads)
p_combinations=0
groups = np.array_split(groupsGlobal[b], no_threads)  # int(no_batches/no_threads)
groupsGlobal=0
pool.map(main, groups)










