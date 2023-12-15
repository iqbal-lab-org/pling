from pathlib import Path
import os
import pandas as pd

def get_pling_root_dir() -> Path:
    return Path(__file__).parent.parent

def get_number_of_batches(outputpath):
    with open(f"{outputpath}/batches/batching_info.txt") as f:
        batch_size, number_of_batches = f.readlines()
    return int(number_of_batches)

def read_in_batch_pairs(filepath):
    pairs=[]
    with open(filepath,"r") as f:
        for line in f:
            genome1, genome2 = (line.strip("[]\n").split(","))
            pairs.append([genome1[1:-1], genome2[2:-1]])
    return pairs

def get_fasta_file_info(genomes_list):
    FASTAFILES_LIST = [el[0] for el in pd.read_csv(genomes_list, header=None).values]
    FASTAFILES = {os.path.splitext(os.path.basename(el))[0]:el for el in FASTAFILES_LIST}
    FASTAEXT = {os.path.splitext(os.path.basename(el))[0]:os.path.splitext(os.path.basename(el))[1] for el in FASTAFILES_LIST}
    FASTAPATH = os.path.dirname(FASTAFILES_LIST[0])
    return FASTAFILES, FASTAEXT, FASTAPATH
