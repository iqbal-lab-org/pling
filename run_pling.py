#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 14:34:47 2023

@author: daria
"""

import subprocess
import argparse
import shutil
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("genomes_list", help="path to list of fasta file paths")
parser.add_argument("output_dir", help="path to output directory")
parser.add_argument("integerisation", choices=["anno", "align"], help="integerisation method: \"anno\" for annotation and \"align\" for alignment")
parser.add_argument("-b","--bakta_db", help="path to bakta database, planned default is where bakta downloads it by default")
parser.add_argument("-j", "--jaccard", default=0.4, help="threshold for initial jaccard network")
parser.add_argument("-d", "--dcj", default=4, help="threshold for final DCJ-indel network")
parser.add_argument("--dedup", default=98.5, help="threshold for separating paralogs in deduplication step (for integerisation from annotation)")
parser.add_argument("--identity", default=80, help="threshold for % of shared sequence between blocks (for integerisation from alignment and for jaccard calculation)")
parser.add_argument("--min_indel_size", default=200, help="minimum size for an indel to be treated as a block (for integerisation from alignment)")
parser.add_argument("--bh_connectivity", default=10, help="minimum number of connections a plasmid need to be considered a blackhole plasmid")
parser.add_argument("--bh_neighbours_edge_density", default=0.2, help="maximum number of edge density between blackhole plasmid neighbours to label the plasmid as blackhole")
parser.add_argument("--small_subcommunity_size_threshold", default=4, help="communities with size up to this parameter will be joined to neighbouring larger subcommunities")
parser.add_argument("--cores", default=1, help="total number of cores/threads")
parser.add_argument("--storetmp", action="store_true", help="don't delete intermediate temporary files")
parser.add_argument("--forceall", action="store_true", help="force snakemake to rerun everything")
parser.add_argument("--timelimit", help="time limit on gurobi")
parser.add_argument("--resources", help="tsv stating number of threads and memory to use for each rule")
parser.add_argument("--profile", help="to run on cluster with corresponding snakemake profile")

args = parser.parse_args()

plingpath = os.path.realpath(os.path.dirname(__file__))

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
    resources = pd.read_csv(f"{plingpath}/resources.tsv", sep="\t")

if not os.path.isdir(f"{args.output_dir}/tmp_files"):
    os.mkdir(f"{args.output_dir}/tmp_files")

configfile = f"{args.output_dir}/tmp_files/config.yaml"
with open(configfile, "w") as config:
    config.write(f"genomes_list: {args.genomes_list}\n\n")
    config.write(f"output_dir: {args.output_dir}\n\n")
    config.write(f"integerisation: {args.integerisation}\n\n")
    config.write(f"bakta_db: {args.bakta_db}\n\n")
    config.write(f"seq_jaccard_threshold: {args.jaccard}\n\n")
    config.write(f"dcj_dist_threshold: {args.dcj}\n\n")
    config.write("prefix: all_plasmids\n\n")
    config.write(f"communities: {args.output_dir}/jaccard/jaccard_communities.txt\n\n")
    config.write(f"dedup_threshold: {args.dedup}\n\n")
    config.write(f"identity_threshold: {args.identity}\n\n")
    config.write(f"length_threshold: {args.min_indel_size}\n\n")
    config.write(f"bh_connectivity: {args.bh_connectivity}\n\n")
    config.write(f"bh_neighbours_edge_density: {args.bh_neighbours_edge_density}\n\n")
    config.write(f"small_subcommunity_size_threshold: {args.small_subcommunity_size_threshold}\n\n")
    config.write(f"timelimit: {timelimit}\n\n")
    for row in resources.index:
        rule = resources.loc[row, "Rule"]
        threads = resources.loc[row, "Threads"]
        mem =  resources.loc[row, "Mem"]
        config.write(f"{rule}_threads: {threads}\n\n")
        config.write(f"{rule}_mem: {mem}\n\n")

#integerisation
if args.integerisation == "anno":
    try:
        print("Building jaccard network...\n")
        subprocess.run(f"snakemake --snakefile {plingpath}/jac_network_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {profile} {forceall}", shell=True, check=True, capture_output=True)
        print("Completed jaccard network.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    try:
        print("Annotating and integerising...\n")
        subprocess.run(f"snakemake --snakefile {plingpath}/anno_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {profile} {forceall}", shell=True, check=True, capture_output=True)
        print("Completed integerisation.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
elif args.integerisation == "align":
    try:
        print("Aligning, intgerising, and building jaccard network...\n")
        subprocess.run(f"snakemake --snakefile {plingpath}/align_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {profile} {forceall}", shell=True, check=True, capture_output=True)
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
    subprocess.run(f"snakemake --snakefile {plingpath}/dcj_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {profile} {forceall}", shell=True, check=True, capture_output=True)
    print("Completed distance calculations and clustering.\n")
except subprocess.CalledProcessError as e:
    print(e.stderr.decode())
    print(e)
    raise e


#delete intermediary files
if args.storetmp == False:
    shutil.rmtree(f"{args.output_dir}/tmp_files")
