import argparse
import os
import math
import numpy as np
import subprocess
from pathlib import Path
from pling.utils import get_fasta_file_info
import itertools

def get_labels(filepath):
    fastafiles, fastaext, fastapath = get_fasta_file_info(filepath)
    genomes = list(fastaext.keys())
    return genomes

def append_pair(smash, smash_threshold, smash_matrix, i, j):
    if not smash:
        return True
    else:
        return (1-smash_matrix[i][j]<=smash_threshold)

def get_pairs(genomes, batch_size, output_dir, containmentpath, smash, smash_matrix = None, smash_threshold = None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    n = len(genomes)
    iter = 0
    batch = -1
    if smash:
        dir = Path(os.path.dirname(containmentpath))
        dir.mkdir(parents=True, exist_ok=True)
        contain_file = open(containmentpath, "w")
    for i,j in itertools.combinations(range(n), 2):
        append = append_pair(smash, smash_threshold, smash_matrix, i, j)
        if append:
            if iter%batch_size==0:
                if iter!=0:
                    batch_file.close()
                batch = batch+1
                batch_file = open(f"{output_dir}/batch_{batch}.txt","w")
            batch_file.write(str([genomes[i], genomes[j]])+"\n")
            iter = iter+1
        elif smash:
            contain_file.write(f"{genomes[i]}\t{genomes[j]}\t{1-smash_matrix[i][j]}\n")
    if smash:
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
        print(e.stderr.decode())
        print(e)
        raise e
    try:
        subprocess.run(f"sourmash compare {sig_path} --max-containment -o {matrixpath}", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def main():
    parser = argparse.ArgumentParser(description="Generate lists of genome pairs per batch")

    parser.add_argument("--genomes_list")
    parser.add_argument("--batch_size", type=int)
    parser.add_argument("--outputpath")
    parser.add_argument("--sourmash", action="store_true")
    parser.add_argument("--smash_threshold",type=float)
    parser.add_argument("--containmentpath")
    parser.add_argument("--dcj_path")

    args = parser.parse_args()

    if args.sourmash:
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

    len_pairs = get_pairs(genomes, args.batch_size, f"{args.outputpath}/batches", args.containmentpath, args.sourmash, smash_matrix, args.smash_threshold)
    number_of_batches = math.ceil(len_pairs/args.batch_size)

    with open(f"{args.outputpath}/batches/batching_info.txt", "w") as f:
        f.write(f"{args.batch_size}\n{number_of_batches}")

if __name__ == "__main__":
    main()
