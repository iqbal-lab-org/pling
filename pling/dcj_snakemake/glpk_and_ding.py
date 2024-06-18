import subprocess
import argparse
from pling.utils import read_in_batch_pairs
from pathlib import Path
from shared_functions import *

def unimog_to_ilp(unimog, lp, genome1, genome2):
    try:
        print("Converting unimog to ILP...\n")
        subprocess.run(f"dingII generate {unimog} -mm --writeilp {lp} -p {genome1} {genome2}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def ilp_GLPK(lp, solution, snakefile_dir, timelimit):
    try:
        print("Solving ILP with GLPK...\n")
        subprocess.run(f"glpsol --lp {lp} -o {solution}.tmp {timelimit}", shell=True, check=True, capture_output=True)
        print("Completed.\n Converting GLPK output to Gurobi output...\n")
        subprocess.run(f"python {snakefile_dir}/glpk_sol_to_gurobi_sol.py <{solution}.tmp >{solution}", shell=True, check=True, capture_output=True)
        print("Completed.\n")
        subprocess.run(f"rm {solution}.tmp", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

def dcj_dist(unimog, solution, genome1, genome2):
    try:
        print("Parsing Gurobi output to get DCJ distances...\n")
        dcj_out = subprocess.run(f"dingII parsesol {unimog} --solgur {solution} -p {genome1} {genome2}", shell=True, check=True, capture_output=True, text=True).stdout
        dist = int(dcj_out.strip().split(" ")[2])
        print("Completed.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    return dist

def batchwise_ding(pairs, containment_distance, containments, integerisation, outputpath, batch, timelimit, snakefile_dir, plasmid_to_community, unimog_path=None):
    dists = []
    for pair in pairs:
        genome1 = pair[0]
        genome2 = pair[1]
        lp = f"ding/ilp/{genome1}~{genome2}.lp"
        solution = f"ding/solutions/{genome1}~{genome2}.sol"
        entry1, entry2 = get_entries(integerisation, genome1, genome2)
        if containments[(genome1,genome2)]<=containment_distance:
            unimog = get_unimog(outputpath, integerisation, plasmid_to_community, batch, genome1, genome2, unimog_path)
            unimog_to_ilp(unimog, lp, entry1, entry2)
            ilp_GLPK(lp, solution, snakefile_dir, timelimit)
            dist = dcj_dist(unimog, solution, entry1, entry2)
            dists.append(f"{genome1}\t{genome2}\t{dist}\n")
    with open(f"{outputpath}/tmp_files/dists_batchwise/batch_{batch}_dcj.tsv", "w") as f:
        for line in dists:
            f.write(line)

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, containment and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--containment_tsv", required=True)
    parser.add_argument("--containment_distance", required=True)
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--communitypath", required=True)
    parser.add_argument("--integerisation", required=True, type=str)
    parser.add_argument("--threads")
    parser.add_argument("--timelimit")
    parser.add_argument("--snakefile_dir", required=True)
    parser.add_argument("--unimog")

    # Parse the arguments
    args = parser.parse_args()

    if args.timelimit == None:
        timelimit=""
    else:
        timelimit=f"--tmlim {args.timelimit}"

    pairs=read_in_batch_pairs(f"{args.outputpath}/batches/batch_{args.batch}.txt")
    containments=get_containment_distances_for_batch(args.containment_tsv)

    output_dirs = [Path(f"ding/ilp"), Path(f"ding/solutions")]
    for dir in output_dirs:
        dir.mkdir(parents=True, exist_ok=True)

    if args.integerisation=="anno":
        plasmid_to_community = get_plasmid_to_community(args.communitypath)
    else:
        plasmid_to_community=None

    if args.integerisation=="skip":
        batchwise_ding(pairs, float(args.containment_distance), containments, args.integerisation, args.outputpath, args.batch, timelimit, args.snakefile_dir, plasmid_to_community, args.unimog)
    else:
        batchwise_ding(pairs, float(args.containment_distance), containments, args.integerisation, args.outputpath, args.batch, timelimit, args.snakefile_dir, plasmid_to_community)

if __name__ == "__main__":
    main()
