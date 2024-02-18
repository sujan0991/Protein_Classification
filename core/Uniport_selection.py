from operator import index
import pandas as pd
import pymol
import pymol.cmd as cmd
import multiprocessing as mp
import numpy as np


def uniport_selection(chain, resi_start, resi_end):

    resi = (resi_start, resi_end)

    s = ""
    if chain and chain != "":

        s = f"chain {chain}"

    if resi and resi is not None and resi != ("", ""):

        s += f" and resi {int(resi[0])}-{int(resi[1])}"

    return s


if __name__ == "__main__":

    ############### create new csv with id and selection from uniprot.csv      ###############

    df = pd.read_csv("uniprot.csv")

    _uniport_selection = np.vectorize(uniport_selection)

    sels = df[["chain", "resi_start", "resi_end"]].fillna("")

    ids = df["pdb_id"].to_numpy()

    selections = _uniport_selection(sels["chain"], sels["resi_start"], sels["resi_end"])

    print("id and selection len", len(ids), len(selections))

    if len(ids) == len(selections):

        d = {"id": ids, "selection": selections}

        df = pd.DataFrame.from_dict(d)

        df.to_csv("uniport_selection.csv", index_label="index")


################
