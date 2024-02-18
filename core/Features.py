## This script will read pdb selection information from a csv file, then extract features
# from PDB files using Pymol

import pymol
import pymol.cmd as cmd
# import psico.editing
import pandas as pd
import sys
from pymol import stored
import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np
import plotly.express as px
import plotly as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import Isomap


start = time.time()
errors_ids = []
diff_aa = []




def get_features(pdb_info):
    
    path = "/Users/Sujan/Desktop/Bio/PDBs/"
    #path = "/group/bioinf_protstr/Ballal/PDBs/"

    # extract PDB id from pdb_info(1st 4 letters)
    sliced = (pdb_info)[:4].lower()

    pdb_name = "pdb" + sliced + ".ent.gz"
    full_path = path + pdb_name
    beforeChain = ""

    if os.path.exists(full_path):
        
        print("str(pdb_info.......", pdb_info)


        try:

            # to save fratures(resi,ss,distance) of the protein
            features_per_pdb = dict()

            # load the file in pymol

            cmd.load(full_path)

            ######     extract selected protein from PDB          ######

            # give a unique name to every protein
            selectedName = "sele_" + sliced

            select_part = ""

            # Some PDB don't have selected part in uniport csv,if that is the case then pdb_info only contains 4 letters
            # check if pdb_info contains selection part
            if len(pdb_info) > 4:

                # check if cath id or pdb id

                beforeChain = pdb_info.split("chain")[0]
                # no cath id without selection part
                if len(beforeChain) > 9:

                    # select_part : selection string 9 letters(after this "1R44 and " part)- cath
                    
                    select_part = (pdb_info)[12:]
                    
                else:

                    # select_part : selection string 9 letters(after this "1R44 and " part)- uniport
                    select_part = (pdb_info)[9:]
                    
                # add Alpha carbon in selection
                selection = select_part + " and n. CA"
                # creaate new element from selected part
                cmd.create(selectedName, selection)

                
                # Delete old PDB, as we already created new element from it
                pdb_name1 = "pdb" + sliced
                cmd.delete(pdb_name1)

            else:

                # if no selection, then only select Alpha carbon
                selection = "n. CA"
                # creaate new element from selected part
                cmd.create(selectedName, selection)

              
                # Delete old PDB, as we already created new element from it
                pdb_name1 = "pdb" + sliced
                cmd.delete(pdb_name1)

            ############.................#############
            
            
            

            # 1.Get_area
            domain_area = cmd.get_area(selectedName)
            print("domain area......", domain_area)
            
            # 2. volume
            
             # to store SS sequence
            sss = []
            sec_struc = []

            # iterate to get ss sequence and store them in sec_struc
            cmd.iterate(selectedName, "sec_struc.append(ss)", space=locals())

            # check if the sequence have any "", if so replace it with 'L' and store the whole sequence in sss
            for c in sec_struc:

                if c.strip() == "":
                    sss.append("L")
                else:
                    sss.append(c)

             
            # number of helices/total number of SS
            
            helices_count = sss.count('H')
            h_volume = helices_count/len(sss)
           
            #3- number of sheets/total number of SS
            sheets_count = sss.count('S')
            s_volume = sheets_count/len(sss)
            
            #4. number of loops/total number of SS
            loops_count = sss.count('L')
            l_volume = loops_count/len(sss)
             
            print("ss........len",selectedName, len(sss),helices_count,sheets_count,loops_count)
            
            
            # 5- number of AAs
            # to store resi sequence
            seq_list = []

            # iterate to get resi sequence and store them in seq_list
            cmd.iterate(selectedName, "seq_list.append([resn])", space=locals())

            # convert list of list to single list
            seq_list = flatten(seq_list)
            print("aa....len", len(seq_list))
            
            #6- number of aromatic residues( phenylalanine, tyrosine, tryptophan and histidine)/number of AAs
            phe_count = seq_list.count('PHE')
            trp_count = seq_list.count('TRP')
            tyr_count = seq_list.count('TYR')
            his_count = seq_list.count('HIS')
            
            total_aromatic_res = phe_count + trp_count + tyr_count + his_count
            
            aromatic_percent = total_aromatic_res/len(seq_list)
            
            print(">>>>>>>total_aromatic_res",total_aromatic_res,aromatic_percent)
            
            #7- number of disulfide linkages (percentage) #1eej
            
            disu_sele = "disulfides_" + sliced
            cmd.select(disu_sele, "CYS/SG and bound_to CYS/SG")
            disulfide_atom_count = cmd.count_atoms(disu_sele)
            disulfide_atom_count_percentage = disulfide_atom_count/len(seq_list)
            print("......disulfide_atom_count_percentage....", disulfide_atom_count_percentage)
            
            # delete selection
            cmd.delete(disu_sele)
            
           
            
            # 10- center of the mass
            com = cmd.centerofmass(selectedName)
            
            
            
            features_per_domain = {}
            cath_id = beforeChain[:7]

            features_per_domain["id"] = cath_id
            features_per_domain["area"] = domain_area
            features_per_domain["h_volume"] = h_volume
            features_per_domain["s_volume"] = s_volume
            features_per_domain["l_volume"] = l_volume
            features_per_domain["aa_count"] = len(seq_list)
            features_per_domain["aromatic_residues_percentage"] = aromatic_percent
            features_per_domain["disulfide_atom_count_percentage"] = disulfide_atom_count_percentage
            features_per_domain["center_of_mass"] = com
        
          
            return features_per_domain
            
            
        except:

            print("error happend")
            # if any error happend, store the PDB id in errors_ids
            errors_ids.append(sliced)

        finally:

            print("finally....")
             # delete the selection, so the new selection may not cause any problem
            cmd.delete(selectedName)
    
    
# conver list of lists to single list
def flatten(t):
    return [item for sublist in t for item in sublist]


# this function return the selection string with PBD id ( pdb_info parameter in read_PDB function)
def get_selection(id, selection):

    s = f"{id}"

    if selection and selection != "":

        s += f" and {selection}"

    return s


if __name__ == "__main__":

    # name_list = ['1R44 and chain A and resi 1-202','1KMD and chain A and resi 8-124','6WLZ and chain C and resi 1-617']

    # read the csv file containing PDB informations
    #df = pd.read_csv("uniport_selection.csv")
    # df = pd.read_csv("cath-selection_class-v4_3_0.csv")

    # # nympy vectorize, to make it faster
    # _get_selection = np.vectorize(get_selection)

    # # select the ids and selection part from df
    # sels = df[["id", "selection"]].fillna("")

    # # call _get_selection to get the selection string with PBD id  and store them in selections
    # selections = _get_selection(sels["id"], sels["selection"])

    # print("selections........", selections[60:70])

    # # sel10k = selections[:10000]

    # # print(".......................",len(sel10k))

    # # to store PDB name (key) with features (value)
    # all_proteins_dictionary = {}

    # # # multiprocessing.Pool()
   
    # with mp.Pool(mp.cpu_count()) as pool:

    #     result = pool.map(get_features, selections)    

    print("selections.........len", len(df.index), len(selections), result)
    
    # output = pd.DataFrame()
    # for (j, value) in enumerate(result):

    #    if value != None:

    #      # print("result..........i",j,value)
          
    #       df_dictionary = pd.DataFrame([value])
    #       output = pd.concat([output, df_dictionary], ignore_index=True)
          
    # print(type(output))
   
    # output.to_pickle("features.pickel")
    
    
    
    
    
    
    
    # creat features_class.csv############
    
    # with open("features.pickel", "rb") as handle:

    #     output = pickle.load(handle)
        
    # df2 = pd.read_csv("cath-domain-description-v4_3_0.csv")
  
    # result = pd.merge(output, df2, on="id")  
    
    # result = result.drop(['index'], axis=1)  
    #result.to_csv("features_class.csv", index_label="index")
    
    ###########