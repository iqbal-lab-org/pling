import argparse
from pling.utils import read_in_batch_pairs
from pathlib import Path
from glpk_and_ding import *

def import_module(module, name):
    import importlib
    try:
        module = importlib.import_module(module)
        return module if name is None else getattr(module, name)
    except ImportError as e:
        msg = f"Missing optional dependency gurobi!"
        raise ValueError(msg) from e

def get_plasmid_to_community(communitypath):
    plasmid_to_community = {}
    with open(communitypath) as communities_fh:
        for community_index, line in enumerate(communities_fh):
            plasmids = line.strip().split()
            for plasmid in plasmids:
                plasmid_to_community[plasmid] = community_index
    return plasmid_to_community

def get_containment_distances_for_batch(containment_tsv):
    containments = {}
    with open(containment_tsv, "r") as f:
        for line in f:
            plasmid_1, plasmid_2, containment = line.strip().split("\t")
            containment = float(containment)
            containments[(plasmid_1,plasmid_2)] = containment
    return containments

def get_unimog(outputpath, integerisation, plasmid_to_community, batch, genome1, genome2, unimog_path=None):
    if integerisation == "align" and unimog_path==None:
        unimog = f"{outputpath}/unimogs/batch_{batch}_align.unimog"
    elif integerisation == "skip":
        unimog=unimog_path
    elif integerisation == "align" and unimog_path!=None:
        unimog=unimog_path
    return unimog

def get_entries(integerisation, genome_1, genome_2):
    if integerisation == "align":
        return f"{genome_1}~{genome_2}:{genome_1}", f"{genome_1}~{genome_2}:{genome_2}"
    elif integerisation == "skip":
        return genome_1, genome_2

def gurobi_flow(pairs, containment_distance, containments, integerisation, outputpath, batch, timelimit, threads, plasmid_to_community, unimog_path):
    compute_DCJ = import_module("gurobi_and_ding", "compute_DCJ")
    dists = []
    for pair in pairs:
        genome1 = pair[0]
        genome2 = pair[1]
        entry1, entry2 = get_entries(integerisation, genome1, genome2)
        if containments[(genome1,genome2)]<=containment_distance:
            unimog = get_unimog(outputpath, integerisation, plasmid_to_community, batch, genome1, genome2, unimog_path)
            dist = compute_DCJ(unimog, entry1, entry2, timelimit, threads)
            dists.append(f"{genome1}\t{genome2}\t{dist}\n")
    return dists

def GLPK_flow(pairs, containment_distance, containments, integerisation, outputpath, batch, timelimit, snakefile_dir, plasmid_to_community, unimog_path):
    output_dirs = [Path(f"ding/ilp"), Path(f"ding/solutions")]
    for dir in output_dirs:
        dir.mkdir(parents=True, exist_ok=True)
    if timelimit == None:
        timelimit=""
    else:
        timelimit=f"--tmlim {args.timelimit}"

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
    return dists

def batchwise_ding(pairs, containment_distance, containments, integerisation, outputpath, distpath, batch, timelimit, threads, snakefile_dir, plasmid_to_community, unimog_path, ilp_solver):
    if ilp_solver == "gurobi":
        dists=gurobi_flow(pairs, containment_distance, containments, integerisation, outputpath, batch, timelimit, threads, plasmid_to_community, unimog_path)
    elif ilp_solver == "GLPK":
        dists=GLPK_flow(pairs, containment_distance, containments, integerisation, outputpath, batch, timelimit, snakefile_dir, plasmid_to_community, unimog_path)
    else:
        raise Exception(f"{ilp_solver} is an invalid solver!")
    with open(distpath, "w") as f:
        for line in dists:
            f.write(line)

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process a pair of genomes and create unimogs, containment and sequence blocks output.")

    # Add the arguments
    parser.add_argument("--batch", required=True, help="Batch number")
    parser.add_argument("--batch_file", required=True, help="Batch list path")
    parser.add_argument("--containment_tsv", required=True)
    parser.add_argument("--containment_distance", required=True)
    parser.add_argument("--outputpath", required=True, help="Path for general output directory")
    parser.add_argument("--distpath", required=True, help="Path to final outputted distance file")
    parser.add_argument("--communitypath", required=True)
    parser.add_argument("--integerisation", required=True, type=str)
    parser.add_argument("--threads")
    parser.add_argument("--timelimit")
    parser.add_argument("--snakefile_dir")
    parser.add_argument("--unimog")
    parser.add_argument("--ilp_solver", required=True)

    # Parse the arguments
    args = parser.parse_args()

    pairs=read_in_batch_pairs(args.batch_file)
    containments=get_containment_distances_for_batch(args.containment_tsv)

    plasmid_to_community=None

    batchwise_ding(pairs, float(args.containment_distance), containments, args.integerisation, args.outputpath, args.distpath, args.batch, args.timelimit, args.threads, args.snakefile_dir, plasmid_to_community, args.unimog, args.ilp_solver)

if __name__ == "__main__":
    main()
