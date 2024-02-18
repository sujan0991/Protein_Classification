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

start = time.time()
errors_ids = []
diff_aa = []


def read_PDB(pdb_info):

    # pdb_info: contain PDB name and selection

    # Path to the PDB files
    path = "/group/bioinf_protstr/Ballal/Unziped_PDB_all/"

    #path = "/Users/Sujan/Desktop/Bio/R_P/"

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
            print('seq_list_fasta...............',selectedName,len(seq_list))
            



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

            print("ss........len", pdb_name, selectedName, len(sss))

            # # save all features in features_per_pdb dict
            features_per_pdb["resns"] = seq_list
            features_per_pdb["ss_struc"] = sss
           
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

            # # if all feature array have same length, return the features_dict
            # if len(seq_list) > 5 and len(seq_list) < 250:

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
    with open("cath_all_proteins_dictionary_SS_AA_with_BM.pickle", "wb") as handle:
        pickle.dump(all_proteins_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        
    
    # ######### convert to one hot vector #########

    
    # main 3 ss
    keys_ss_struc = ["H", "S", "L"]
    # to store features with PDB name in one hot format
    final_features = dict()

    for k in all_proteins_dictionary.keys():

        # GET THE values from all_proteins_dictionary
        val = all_proteins_dictionary[k]
        # to store features in one hot format
        features = []

        for k2 in range(len(val["ss_struc"])):

            featuresDic_ss_struc = dict.fromkeys(
                keys_ss_struc, 0
            )  # dictionary of ss name, the form of one hot vector. First, set all of them  to 0

                # if the ss_struc matches the name in keys_ss_struc, set that bit to 1, rest will be 0
            if val["ss_struc"][k2] in keys_ss_struc:
                    featuresDic_ss_struc[val["ss_struc"][k2]] = 1

               
                # append all three features in features
            features.append(list(featuresDic_ss_struc.values()))

            
        # store all feature in final_features where the key is the PDB id and value id features list
        final_features[k] = features

    print("final_features.....len", len(final_features))
    print("have different aa", len(set(diff_aa)), set(diff_aa))
    print("errors_ids....", errors_ids)
    
    # for key, value in final_features.items():

    #     print(">>>>>>>>>>>>","key:", key,"value len:", len(value))

    # save final_features in pickle
    with open("cath_final_features_with_BM_only_SS.pickle", "wb") as handle:
        pickle.dump(final_features, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(time.time() - start)





