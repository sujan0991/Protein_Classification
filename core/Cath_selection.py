from operator import index
import pandas as pd
import pymol
import pymol.cmd as cmd
import multiprocessing as mp
import numpy as np


def cath_selection(id, selection):

    if len(id) == 7:

        pdb_id = id[:4]
        chain = id[4]
        domain = id[5:]

        #s = f"{pdb_id}"
        if chain and chain != "":

            s = f"chain {chain}"

        if selection and selection is not None:

            s += f" and resi {selection}"

        return s


if __name__ == "__main__":

    # ########### convert cath txt file to csv    ################

    df = pd.read_csv(
        "cath-domain-boundaries-seqreschopping-v4_3_0.txt", delim_whitespace=True
    )

    df.to_csv(
        "cath-domain-boundaries-seqreschopping-v4_3_0.csv",
        index=True,
        header=["id", "selection"],
        index_label="index",
    )

    #############

    ####### create new cath csv with id and selection   ##########

    df = pd.read_csv("cath-domain-boundaries-seqreschopping-v4_3_0.csv")
    _cath_selection = np.vectorize(cath_selection)
    sels = df[["id", "selection"]].fillna("")
    ids = df["id"].to_numpy()
    pdb_ids = []

    for single_id in ids:
        pdb_id = single_id[:4]
        pdb_ids.append(pdb_id)
        
    print("61------67",ids[60:70], pdb_ids[60:70])  
    selections = _cath_selection(sels["id"], sels["selection"])
    print("id and selection len", len(ids),len(pdb_ids), len(selections))

    if len(ids) == len(selections):
        d = {"id": ids,"pdb_id": pdb_ids, "selection": selections}
        df = pd.DataFrame.from_dict(d)
        df.to_csv("cath-selection-v4_3_0.csv", index_label="index")

    # #############


# marge cath-selection-v4_3_0.csv and cath-domain-description-v4_3_0.csv
#   df=pd.read_csv("cath-selection-v4_3_0.csv")
#   df2 = pd.read_csv("cath-domain-description-v4_3_0.csv")
  
#   result = pd.merge(df, df2, on="id")
#   result.to_csv("cath-selection_class-v4_3_0.csv", index_label="index")

  
#   print("???????????",result[:10])