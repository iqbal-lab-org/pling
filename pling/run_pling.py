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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes_list", help="path to list of fasta file paths")
    parser.add_argument("output_dir", help="path to output directory")
    parser.add_argument("integerisation", choices=["anno", "align"], help="integerisation method: \"anno\" for annotation and \"align\" for alignment")
    parser.add_argument("-b","--bakta_db", help="path to bakta database, planned default is where bakta downloads it by default")
    parser.add_argument("-j", "--jaccard-distance", default=0.6, help="threshold for initial jaccard network")
    parser.add_argument("-d", "--dcj", default=4, help="threshold for final DCJ-indel network")
    parser.add_argument("--dedup", action="store_true", help="whether or not the deduplicate when integerising from annotation")
    parser.add_argument("--dedup_threshold", default=98.5, help="threshold for separating paralogs in deduplication step (for integerisation from annotation)")
    parser.add_argument("--identity", default=80, help="threshold for percentage of shared sequence between blocks (for integerisation from alignment and for jaccard calculation)")
    parser.add_argument("--min_indel_size", default=200, help="minimum size for an indel to be treated as a block (for integerisation from alignment)")
    parser.add_argument("--bh_connectivity", default=10, help="minimum number of connections a plasmid need to be considered a blackhole plasmid")
    parser.add_argument("--bh_neighbours_edge_density", default=0.2, help="maximum number of edge density between blackhole plasmid neighbours to label the plasmid as blackhole")
    parser.add_argument("--small_subcommunity_size_threshold", default=4, help="communities with size up to this parameter will be joined to neighbouring larger subcommunities")
    parser.add_argument("--cores", default=1, help="total number of cores/threads")
    parser.add_argument("--storetmp", action="store_true", help="don't delete intermediate temporary files")
    parser.add_argument("--forceall", action="store_true", help="force snakemake to rerun everything")
    parser.add_argument("--ilp_solver", choices=["GLPK", "gurobi"], default="GLPK",
                        help="ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. "
                             "If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH.")
    parser.add_argument("--timelimit", help="time limit in seconds for ILP solver")
    parser.add_argument("--resources", help="tsv stating number of threads and memory to use for each rule")
    parser.add_argument("--profile", help="to run on cluster with corresponding snakemake profile")

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
    with open(configfile, "w") as config:
        config.write(f"genomes_list: {args.genomes_list}\n\n")
        config.write(f"output_dir: {args.output_dir}\n\n")
        config.write(f"integerisation: {args.integerisation}\n\n")
        config.write(f"bakta_db: {args.bakta_db}\n\n")
        config.write(f"seq_jaccard_distance: {args.jaccard_distance}\n\n")
        config.write(f"dcj_dist_threshold: {args.dcj}\n\n")
        config.write("prefix: all_plasmids\n\n")
        config.write(f"communities: {args.output_dir}/jaccard/jaccard_communities\n\n")
        if args.dedup:
            config.write(f"dedup: {args.dedup}\n\n")
            config.write(f"dedup_threshold: {args.dedup_threshold}\n\n")
        config.write(f"identity_threshold: {args.identity}\n\n")
        config.write(f"length_threshold: {args.min_indel_size}\n\n")
        config.write(f"bh_connectivity: {args.bh_connectivity}\n\n")
        config.write(f"bh_neighbours_edge_density: {args.bh_neighbours_edge_density}\n\n")
        config.write(f"small_subcommunity_size_threshold: {args.small_subcommunity_size_threshold}\n\n")
        config.write(f"ilp_solver: {args.ilp_solver}\n\n")
        config.write(f"timelimit: {timelimit}\n\n")
        for row in resources.index:
            rule = resources.loc[row, "Rule"]
            threads = resources.loc[row, "Threads"]
            mem =  resources.loc[row, "Mem"]
            config.write(f"{rule}_threads: {threads}\n\n")
            config.write(f"{rule}_mem: {mem}\n\n")

    return configfile, tmp_dir, forceall, profile


def pling(args):
    configfile, tmp_dir, forceall, profile = make_config_file(args)

    snakemake_args = f"--configfile {configfile} --cores {args.cores} --use-conda --rerun-incomplete {profile} {forceall}"

    #integerisation
    if args.integerisation == "anno":
        try:
            print("Building jaccard network...\n")
            subprocess.run(f"snakemake --snakefile {get_pling_path()}/jac_network_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
            print("Completed jaccard network.\n")
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
            print("Aligning, intgerising, and building jaccard network...\n")
            subprocess.run(f"snakemake --snakefile {get_pling_path()}/align_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
            print("Completed integerisation and jaccard network.\n")
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
    if not args.storetmp:
        shutil.rmtree(tmp_dir)


def main():
    args = parse_args()
    pling(args)

if __name__ == "__main__":
    main()
