"""
Created on Tue Aug 15 14:34:47 2023

@author: daria
"""

import subprocess
import argparse
import shutil
import os
import pandas as pd
from pathlib import Path
import subprocess
import yaml

def get_version():
    try:
        pling_dir = get_pling_path()
        output = subprocess.check_output("git describe --tags", shell=True, cwd=pling_dir).strip()
        output = output.decode("utf-8").split('-')[0]
    except:
        output = "Error getting version."
    return output


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version=get_version())
    parser.add_argument("genomes_list", help="Path to list of fasta file paths.")
    parser.add_argument("output_dir", help="Path to output directory.")
    parser.add_argument("integerisation", choices=["anno", "align"],
                        help="Integerisation method: \"anno\" for annotation and \"align\" for alignment."
                             " Note that recommended method is integerisation from alignment.")
    parser.add_argument("--containment_distance", default=0.5, help="Threshold for initial containment network.")
    parser.add_argument("--dcj", default=4, help="Threshold for final DCJ-Indel network.")
    parser.add_argument("--batch_size", default = 200, help="How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps).")
    parser.add_argument("--sourmash", action="store_true", help="Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.")
    parser.add_argument("--sourmash_threshold", default=0.85, help="Threshold for filtering with sourmash.")
    parser.add_argument("--identity", default=80, help="Threshold for percentage of shared sequence between blocks (for integerisation from alignment and for containment calculation).")
    parser.add_argument("--min_indel_size", default=200, help="Minimum size for an indel to be treated as a block (for integerisation from alignment).")
    parser.add_argument("--bh_connectivity", default=10, help="Minimum number of connections a plasmid need to be considered a hub plasmid.")
    parser.add_argument("--bh_neighbours_edge_density", default=0.2, help="Maximum number of edge density between hub plasmid neighbours to label the plasmid as hub.")
    parser.add_argument("--small_subcommunity_size_threshold", default=4, help="Communities with size up to this parameter will be joined to neighbouring larger subcommunities.")
    parser.add_argument("--plasmid_metadata", help="Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list.")
    parser.add_argument("--ilp_solver", choices=["GLPK", "gurobi"], default="GLPK",
                        help="ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. "
                             "If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH.")
    parser.add_argument("--timelimit", help="Time limit in seconds for ILP solver.")
    parser.add_argument("--resources", help="tsv stating number of threads and memory to use for each rule.")
    parser.add_argument("--cores", default=1, help="Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.)")
    parser.add_argument("--profile", help="To run on a cluster with corresponding snakemake profile.")
    #parser.add_argument("--storetmp", action="store_true", help="Don't delete intermediate temporary files.")
    parser.add_argument("--forceall", action="store_true", help="Force snakemake to rerun everything.")
    parser.add_argument("--dedup", action="store_true", help="Whether or not to deduplicate (for integerisation from annotation).")
    parser.add_argument("--dedup_threshold", default=98.5, help="Threshold for separating paralogs in deduplication step (for integerisation from annotation).")
    parser.add_argument("--bakta_db", help="Path to bakta database (required for integerisation from annotation).")

    args = parser.parse_args()
    return args


def get_pling_path():
    plingpath = os.path.realpath(os.path.dirname(__file__))
    return plingpath

def make_config_file(args):
    #make configfile
    forceall = ""
    if args.forceall == True:
        forceall = "--forceall"

    if args.timelimit==None:
        timelimit = "None"
    else:
        timelimit= args.timelimit

    if args.plasmid_metadata==None:
        metadata = "None"
    else:
        metadata= args.metadata

    profile = ""
    if args.profile!=None:
        profile = f"--profile {args.profile}"
    if args.resources!=None:
        resources = pd.read_csv(args.resources, sep="\t")
    else:
        resources = pd.read_csv(f"{get_pling_path()}/resources.tsv", sep="\t")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = output_dir/"tmp_files"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    configfile = f"{args.output_dir}/tmp_files/config.yaml"
    config_dict = {"genomes_list": str(args.genomes_list), "output_dir": str(args.output_dir), "integerisation": str(args.integerisation), "bakta_db": str(args.bakta_db), "seq_containment_distance": float(args.containment_distance), "dcj_dist_threshold": int(args.dcj), "prefix": "all_plasmids","communities": f"{args.output_dir}/containment/containment_communities", "identity_threshold": float(args.identity), "length_threshold": int(args.min_indel_size), "bh_connectivity": int(args.bh_connectivity), "bh_neighbours_edge_density": float(args.bh_neighbours_edge_density), "small_subcommunity_size_threshold": int(args.small_subcommunity_size_threshold),"metadata": metadata, "ilp_solver": str(args.ilp_solver), "timelimit": timelimit, "batch_size": int(args.batch_size)}
    if args.dedup:
        config_dict["dedup"] = str(args.dedup)
        config_dict["dedup_threshold"] = str(args.dedup_threshold)
    if args.sourmash:
        config_dict["sourmash"] = str(args.sourmash)
        config_dict["sourmash_threshold"] = str(args.sourmash_threshold)
    for row in resources.index:
        rule = resources.loc[row, "Rule"]
        threads = resources.loc[row, "Threads"]
        mem =  resources.loc[row, "Mem"]
        if not pd.isna(threads):
            config_dict[f"{rule}_threads"] = int(threads)
        config_dict[f"{rule}_mem"] = int(mem)

    with open(configfile, 'w') as config:
        yaml.dump(config_dict, config)

    return configfile, tmp_dir, forceall, profile

def pling(args):
    configfile, tmp_dir, forceall, profile = make_config_file(args)

    snakemake_args = f"--configfile {configfile} --cores {args.cores} --use-conda --rerun-incomplete --nolock {profile} {forceall}"

    #batching
    try:
        print("Batching...\n")
        subprocess.run(f"snakemake --snakefile {get_pling_path()}/batching/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
        print("Completed batching.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

    #integerisation
    if args.integerisation == "anno":
        try:
            print("Building containment network...\n")
            subprocess.run(f"snakemake --snakefile {get_pling_path()}/jac_network_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
            print("Completed containment network.\n")
        except subprocess.CalledProcessError as e:
            print(e.stderr.decode())
            print(e)
            raise e
        try:
            print("Annotating and integerising...\n")
            subprocess.run(f"snakemake --snakefile {get_pling_path()}/anno_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
            print("Completed integerisation.\n")
        except subprocess.CalledProcessError as e:
            print(e.stderr.decode())
            print(e)
            raise e
    elif args.integerisation == "align":
        try:
            print("Aligning, integerising, and building containment network...\n")
            subprocess.run(f"snakemake --snakefile {get_pling_path()}/align_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
            print("Completed integerisation and containment network.\n")
        except subprocess.CalledProcessError as e:
            print(e.stderr.decode())
            print(e)
            raise e
    else:
        raise Exception(f"\"{args.integerisation}\" is not a valid integerisation method!")

    #ding and clustering
    try:
        print("Calculating DCJ-Indel distance and clustering...\n")
        subprocess.run(f"snakemake --snakefile {get_pling_path()}/dcj_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
        print("Completed distance calculations and clustering.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e


    #delete intermediary files
    shutil.rmtree(tmp_dir)


def main():
    args = parse_args()
    pling(args)

if __name__ == "__main__":
    main()
