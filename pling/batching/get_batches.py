import argparse
import os
import math
import numpy as np
import pandas as pd
import subprocess
from pathlib import Path
from pling.utils import get_fasta_file_info
import itertools
import sys

def get_labels(filepath):
    fastafiles, fastaext, fastapath = get_fasta_file_info(filepath)
    genomes = list(fastaext.keys())
    return genomes

def append_pair(smash, smash_threshold, smash_matrix, i, j):
    if not smash:
        return True
    else:
        return (1-smash_matrix[i][j]<=smash_threshold)
    
def previous_pair(genome_1, genome_2, prev_genomes):
    for genomes in prev_genomes:
        if genome_1 in genomes and genome_2 in genomes:
            return True
    return False

def get_pairs(genomes, batch_size, output_dir, containmentpath, smash, smash_only, smash_matrix = None, smash_threshold = None, prev_genomes = []):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    n = len(genomes)
    iter = 0
    batch = -1
    if smash or smash_only:
        dir = Path(os.path.dirname(containmentpath))
        dir.mkdir(parents=True, exist_ok=True)
        contain_file = open(containmentpath, "w")
        contain_file.write("plasmid_1\tplasmid_2\tdistance\n")
    for i,j in itertools.combinations(range(n), 2):
        if not previous_pair(genomes[i], genomes[j], prev_genomes):
            append = append_pair(smash, smash_threshold, smash_matrix, i, j)
            if append:
                if iter%batch_size==0:
                    if iter!=0:
                        batch_file.close()
                    batch = batch+1
                    batch_file = open(f"{output_dir}/batch_{batch}.txt","w")
                batch_file.write(str([genomes[i], genomes[j]])+"\n")
                iter = iter+1
                if smash_only:
                    contain_file.write(f"{genomes[i]}\t{genomes[j]}\t{1-smash_matrix[i][j]}\n")
            elif smash:
                contain_file.write(f"{genomes[i]}\t{genomes[j]}\t{1-smash_matrix[i][j]}\n")
    if smash or smash_only:
        contain_file.close()
    batch_file.close()
    return iter

def containment_file(not_pairs, genome_index, smash_matrix, containmentpath):
    dir = Path(os.path.dirname(containmentpath))
    dir.mkdir(parents=True, exist_ok=True)
    with open(containmentpath, "w") as f:
        for el in not_pairs:
            i = genome_index[el[0]]
            j = genome_index[el[1]]
            f.write(f"{el[0]}\t{el[1]}\t{1-smash_matrix[i][j]}\n")

def run_smash(genome_list, sig_path, matrixpath):
    try:
        subprocess.run(f"sourmash sketch dna --from-file {genome_list} -o {sig_path}", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode(), file=sys.stderr)
        print(e)
        raise e
    try:
        subprocess.run(f"sourmash compare {sig_path} --max-containment -o {matrixpath}", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode(), file=sys.stderr)
        print(e)
        raise e

def main():
    parser = argparse.ArgumentParser(description="Generate lists of genome pairs per batch")

    parser.add_argument("--genomes_list")
    parser.add_argument("--batch_size", type=int)
    parser.add_argument("--outputpath")
    parser.add_argument("--sourmash", action="store_true")
    parser.add_argument("--sourmash_only", action="store_true")
    parser.add_argument("--smash_threshold",type=float)
    parser.add_argument("--containmentpath")
    parser.add_argument("--prev_typing", nargs="*")

    args = parser.parse_args()

    if args.sourmash or args.sourmash_only:
        sig_dir = Path(f"{args.outputpath}/sourmash")
        sig_dir.mkdir(parents=True, exist_ok=True)
        sig_path = sig_dir/"all_plasmids.sig"
        matrixpath = sig_dir/"smash_containment_matrix"
        run_smash(args.genomes_list, sig_path, matrixpath)
        genomes = get_labels(f"{matrixpath}.labels.txt")
        smash_matrix = np.load(matrixpath, mmap_mode='r')
    else:
        genomes = get_labels(args.genomes_list)
        smash_matrix = None

    if args.prev_typing:
        prev_genomes = [pd.read_csv(f"{path}/typing.tsv", sep="\t")["plasmid"].to_list() for path in args.prev_typing]
        prev_hubs = [pd.read_csv(f"{path}/hub_plasmids.csv", sep="\t")["hub_plasmids"].to_list() for path in args.prev_typing]
        prev_genomes = [prev_genomes[i]+prev_hubs[i] for i in range(len(args.prev_typing))]
    else:
        prev_genomes=[]

    len_pairs = get_pairs(genomes, args.batch_size, f"{args.outputpath}/batches", args.containmentpath, args.sourmash, args.sourmash_only, smash_matrix, args.smash_threshold,prev_genomes=prev_genomes)
    number_of_batches = math.ceil(len_pairs/args.batch_size)

    with open(f"{args.outputpath}/batches/batching_info.txt", "w") as f:
        f.write(f"{args.batch_size}\n{number_of_batches}")

if __name__ == "__main__":
    main()
