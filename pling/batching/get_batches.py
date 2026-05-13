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

def get_labels(filepath, paths=False):
    fastafiles, fastaext, fastapath = get_fasta_file_info(filepath)
    genomes = list(fastaext.keys())
    if paths:
        return genomes, fastafiles
    return genomes

def append_pair(smash, smash_threshold, smash_matrix, i, j):
    if not smash:
        return True
    else:
        return (1-smash_matrix[i][j]<=smash_threshold)
    
def write_content(list_a, list_b, index_map, output_dir, batch_size, smash, smash_only, smash_threshold, smash_matrix, contain_file, batch=-1, iter=0):
    for genome_1, genome_2 in itertools.combinations(list_a, list_b):
        i = index_map[genome_1]
        j = index_map[genome_2]
        append = append_pair(smash, smash_threshold, smash_matrix, i, j)
        if append:
            if iter%batch_size==0:
                if iter!=0:
                    batch_file.close()
                batch = batch+1
                batch_file = open(f"{output_dir}/batch_{batch}.txt","w")
            batch_file.write(str([genome_1, genome_2])+"\n")
            iter = iter+1
            if smash_only:
                contain_file.write(f"{genome_1}\t{genome_2}\t{1-smash_matrix[i][j]}\n")
        elif smash:
            contain_file.write(f"{genome_1}\t{genome_2}\t{1-smash_matrix[i][j]}\n")
    return contain_file, batch_file, batch, iter

def get_pairs(genomes, batch_size, output_dir, containmentpath, smash, smash_only, smash_matrix = None, smash_threshold = None, prev_genomes = [], new_genomes=[]):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    index_map = {genome: i for i, genome in enumerate(genomes)}
    if smash or smash_only:
        dir = Path(os.path.dirname(containmentpath))
        dir.mkdir(parents=True, exist_ok=True)
        contain_file = open(containmentpath, "w")
        contain_file.write("plasmid_1\tplasmid_2\tdistance\n")
    for list_a, list_b in itertools.combinations(prev_genomes):
        contain_file, batch_file, batch, iter = write_content(list_a, list_b, index_map, output_dir, batch_size, smash, smash_only, smash_threshold, smash_matrix, contain_file)
    if prev_genomes == []:
        new_genomes = genomes
    if new_genomes!=[]:
        contain_file, batch_file, batch, iter = write_content(new_genomes, new_genomes, index_map, output_dir, batch_size, smash, smash_only, smash_threshold, smash_matrix, contain_file, batch=batch, iter=iter)
    if smash or smash_only:
        contain_file.close()
    batch_file.close()
    return iter

def run_smash_w_prev(genomes, prev_signatures, fastafiles, sig_path, matrixpath):
    if genomes!=[]:
        try:
            paths = " ".join([fastafiles[genome] for genome in genomes])
            subprocess.run(f"sourmash sketch dna {paths} -o {sig_path}", shell=True, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR IN RULE: SOURMASH SKETCH FAILING WITH NEW GENOMES", file=sys.stderr)
            print(e.stderr.decode(), file=sys.stderr)
            print("END ERROR MSG", file=sys.stderr)
            raise e
    else:
        sig_path = ""
    try:
        prev_sigs=" ".join([sig for sig in prev_signatures])
        subprocess.run(f"sourmash compare {sig_path} {prev_sigs} --max-containment -o {matrixpath}", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR IN RULE: SOURMASH COMPARE FAILING", file=sys.stderr)
        print(e.stderr.decode(), file=sys.stderr)
        print("END ERROR MSG", file=sys.stderr)
        raise e

def run_smash(genome_list, sig_path, matrixpath):
    try:
        subprocess.run(f"sourmash sketch dna --from-file {genome_list} -o {sig_path}", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR IN RULE: SOURMASH SKETCH FAILING WITH {genome_list}", file=sys.stderr)
        print(e.stderr.decode(), file=sys.stderr)
        print("END ERROR MSG", file=sys.stderr)
        raise e
    try:
        subprocess.run(f"sourmash compare {sig_path} --max-containment -o {matrixpath}", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR IN RULE: SOURMASH COMPARE FAILING WITH {genome_list}", file=sys.stderr)
        print(e.stderr.decode(), file=sys.stderr)
        print("END ERROR MSG", file=sys.stderr)
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

    if args.prev_typing:
        genomes, fastafiles = get_labels(args.genomes_list, paths=True)
        prev_genomes = [pd.read_csv(f"{path}/typing.tsv", sep="\t")["plasmid"].to_list() for path in args.prev_typing]
        prev_hubs = [pd.read_csv(f"{path}/hub_plasmids.csv", sep="\t")["hub_plasmids"].to_list() for path in args.prev_typing]
        prev_genomes = [prev_genomes[i]+prev_hubs[i] for i in range(len(args.prev_typing))]
        new_genomes = list(set(genomes)-set(prev_genomes))
        if args.sourmash or args.sourmash_only:
                prev_signatures = [f"{path}/sourmash/all_plasmids.sig" for path in args.prev_typing]
    else:
        prev_genomes=[]

    if args.sourmash or args.sourmash_only:
        sig_dir = Path(f"{args.outputpath}/sourmash")
        sig_dir.mkdir(parents=True, exist_ok=True)
        sig_path = sig_dir/"all_plasmids.sig"
        matrixpath = sig_dir/"smash_containment_matrix"
        if args.prev_typing:
            run_smash(new_genomes, prev_signatures, fastafiles, matrixpath, prev_signatures)
        else:
            run_smash(args.genomes_list, sig_path, matrixpath)
        genomes = get_labels(f"{matrixpath}.labels.txt")
        smash_matrix = np.load(matrixpath, mmap_mode='r')
    else:
        genomes = get_labels(args.genomes_list)
        smash_matrix = None

    len_pairs = get_pairs(genomes, args.batch_size, f"{args.outputpath}/batches", args.containmentpath, args.sourmash, args.sourmash_only, smash_matrix, args.smash_threshold,prev_genomes=prev_genomes)
    number_of_batches = math.ceil(len_pairs/args.batch_size)

    with open(f"{args.outputpath}/batches/batching_info.txt", "w") as f:
        f.write(f"{args.batch_size}\n{number_of_batches}")

if __name__ == "__main__":
    main()
