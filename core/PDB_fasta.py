## This script will read pdb selection information from a csv file, then extract features
# from PDB files using Pymol

import pymol
import pymol.cmd as cmd
import pandas as pd
import sys
from pymol import stored
import multiprocessing as mp
import os.path
import pickle
import time
import numpy as np

start = time.time()
errors_ids = []
diff_aa = []


def read_PDB(pdb_info):

    # pdb_info: contain PDB name and selection

    # Path to the PDB files
    path = "/group/bioinf_protstr/Ballal/PDBs/"

    

    # extract PDB id from pdb_info(1st 4 letters)
    sliced = (pdb_info)[:4].lower()

    pdb_name = "pdb" + sliced + ".ent.gz"
    full_path = path + pdb_name
    beforeChain = ""

    if os.path.exists(full_path):

        print("str(pdb_info.......", full_path, pdb_info)

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

                print("selectedName.....", selectedName)
                print("selection.....", selection)

                # Delete old PDB, as we already created new element from it
                pdb_name1 = "pdb" + sliced
                cmd.delete(pdb_name1)

            else:

                # if no selection, then only select Alpha carbon
                selection = "n. CA"
                # creaate new element from selected part
                cmd.create(selectedName, selection)

                print("selectedName.....", selectedName)
                print("selection.....?????????????????", selection)

                # Delete old PDB, as we already created new element from it
                pdb_name1 = "pdb" + sliced
                cmd.delete(pdb_name1)

            ############.................#############



            # to store resi sequence
            seq_list = []
            
            seq_list = cmd.get_fastastr(selectedName)
            print(pdb_name,"-----------",seq_list)

            features_per_pdb["resns"] = seq_list
            

            # create a new dictionary with pdb name as key
            features_dict = {}
            # temp_dic[sliced] = seq_list.values.tolist()

            # if cath csv, then save with cath id
            if len(beforeChain) > 9:

                cath_id = beforeChain[:7]

                features_dict[cath_id] = features_per_pdb
            # else save with pdb id
            else:

                features_dict[sliced] = features_per_pdb

            # delete the selection, so the new selection may not cause any problem
            cmd.delete(selectedName)

            
                
            if len(seq_list) > 53 and len(seq_list) < 249:

                  return features_dict

        except:

            print("error happend")
            # if any error happend, store the PDB id in errors_ids
            errors_ids.append(sliced)

        finally:

            print("finally....")


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
    df = pd.read_csv("cath-selection-v4_3_0.csv")

    # nympy vectorize, to make it faster
    _get_selection = np.vectorize(get_selection)

    # select the ids and selection part from df
    sels = df[["id", "selection"]].fillna("")

    # call _get_selection to get the selection string with PBD id  and store them in selections
    selections = _get_selection(sels["id"], sels["selection"])

    print("selections........", selections[60:70])

    # sel10k = selections[:10000]

    # print(".......................",len(sel10k))

    # to store PDB name (key) with features (value)
    all_proteins_dictionary = {}

    # # multiprocessing.Pool()
    # run read_PDB in parallel for all PDBS using multiprocessing and store the results in result
    with mp.Pool(mp.cpu_count()) as pool:

        result = pool.map(read_PDB, selections)

    print("selections.........len", len(df.index), len(selections), len(result))

    # loop over the result and store every PDB features in all_proteins_dictionary as value
    for (j, value) in enumerate(result):

        # print("result..........i",j,value)

        if value != None:

            all_proteins_dictionary[list(value.keys())[0]] = list(value.values())[0]

    print("all_proteins_dictionary.....len", len(all_proteins_dictionary))

    # save all_proteins_dictionary in pickle
    with open("cath_all_proteins_dictionary_fasta.pickle", "wb") as handle:
        pickle.dump(all_proteins_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    
  


