import argparse
import os
import math
import numpy as np
import subprocess
from pathlib import Path
from utils import get_fasta_file_info

def write_batch_file(output_dir, batches, batch):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(f"{output_dir}/batch_{batch}.txt","w") as batch_list:
        batch_list.write("\n".join([str(el) for el in batches[str(batch)]]))

def get_labels(filepath):
    fastafiles, fastaext, fastapath = get_fasta_file_info(filepath)
    genomes = list(fastaext.keys())
    genome_index = {genome:i for i, genome in enumerate(genomes)}
    return genomes, genome_index

def append_pair(smash, smash_threshold, smash_matrix, i, j):
    if not smash:
        return True
    else:
        return (1-smash_matrix[i][j]<=smash_threshold)

def get_pairs(genomes, smash, smash_matrix = None, smash_threshold = None):
    genome_pairs=[]
    not_pairs=[]
    n = len(genomes)
    for i in range(n):
        j=0
        while j<i:
            append = append_pair(smash, smash_threshold, smash_matrix, i, j)
            if append:
                genome_pairs.append([genomes[i], genomes[j]])
            else:
                not_pairs.append([genomes[i], genomes[j]])
            j=j+1
    return genome_pairs, not_pairs

def jaccard_file(not_pairs, genome_index, smash_matrix, jaccardpath):
    dir = Path(os.path.dirname(jaccardpath))
    dir.mkdir(parents=True, exist_ok=True)
    with open(jaccardpath, "w+") as f:
        for el in not_pairs:
            i = genome_index[el[0]]
            j = genome_index[el[1]]
            f.write(f"{el[0]}\t{el[1]}\t{1-smash_matrix[i][j]}\n")

def dcj_file(not_pairs, genomes, dcj_path):
    dir = Path(os.path.dirname(dcj_path))
    dir.mkdir(parents=True, exist_ok=True)
    with open(dcj_path, "w+") as f:
        for el in not_pairs:
            f.write(f"{el[0]}\t{el[1]}\n")
        for genome in genomes:
            f.write(f"{genome}\t{genome}\t0\n")

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
    parser.add_argument("--jaccardpath")
    parser.add_argument("--dcj_path")

    args = parser.parse_args()

    if args.sourmash:
        sig_dir = Path(f"{args.outputpath}/sourmash")
        sig_dir.mkdir(parents=True, exist_ok=True)
        sig_path = sig_dir/"all_plasmids.sig"
        matrixpath = sig_dir/"smash_jaccard_matrix"
        run_smash(args.genomes_list, sig_path, matrixpath)
        genomes, genome_index = get_labels(f"{matrixpath}.labels.txt")
        smash_matrix = np.load(matrixpath)
    else:
        genomes, genome_index = get_labels(args.genomes_list)
        smash_matrix = None

    pairs, not_pairs = get_pairs(genomes, args.sourmash, smash_matrix, args.smash_threshold)

    if args.sourmash:
        jaccard_file(not_pairs, genome_index, smash_matrix, args.jaccardpath)

    dcj_file(not_pairs, genomes, args.dcj_path)

    number_of_batches = math.ceil(len(pairs)/args.batch_size)
    batches = {}
    for i in range(number_of_batches-1):
        batches[str(i)] = [pairs[j] for j in range(i*args.batch_size, (i+1)*args.batch_size)]
        write_batch_file(f"{args.outputpath}/batches", batches, str(i))
    batches[str(number_of_batches-1)] = [pairs[j] for j in range((number_of_batches-1)*args.batch_size, len(pairs))]
    write_batch_file(f"{args.outputpath}/batches", batches, str(number_of_batches-1))

    with open(f"{args.outputpath}/batches/batching_info.txt", "w+") as f:
        f.write(f"{args.batch_size}\n{number_of_batches}")

if __name__ == "__main__":
    main()
