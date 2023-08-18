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

parser = argparse.ArgumentParser()
parser.add_argument("genomes_list", help="path to list of fasta file names")
parser.add_argument("output_dir", help="path to output directory")
parser.add_argument("integerisation", choices=["anno", "align"], help="integerisation method: \"anno\" for annotation and \"align\" for alignment")
parser.add_argument("-b","--bakta_db", help="path to bakta database, planned default is where bakta downloads it by default")
parser.add_argument("-j", "--jaccard", default=0.4, help="threshold for initial jaccard network")
parser.add_argument("-d", "--dcj", default=4, help="threshold for final DCJ-indel network")
parser.add_argument("--dedup", default=98.5, help="threshold for separating paralogs in deduplication step")
parser.add_argument("--identity", default=80, help="threshold for % of shared sequence between blocks")
parser.add_argument("--cores", default=1, help="number of cores/threads")
parser.add_argument("--storetmp", action="store_true", help="don't delete intermediate temporary files")
parser.add_argument("--forceall", action="store_true", help="force snakemake to rerun everything")
parser.add_argument("--timelimit", help="time limit on gurobi")

args = parser.parse_args()

plingpath = os.path.realpath(os.path.dirname(__file__))

forceall = ""
if args.forceall == True:
    forceall = "--forceall"
    
if args.timelimit==None:
    timelimit = "None"
else:
    timelimit= args.timelimit

#make configfile
if not os.path.isdir(f"{args.output_dir}/tmp_files"):
    os.mkdir(f"{args.output_dir}/tmp_files")

configfile = f"{args.output_dir}/tmp_files/config.yaml"
with open(configfile, "w") as config:
    config.write(f"genomes_list: {args.genomes_list}\n\n")
    config.write(f"fastas: {args.fasta_dir}\n\n")
    config.write(f"output_dir: {args.output_dir}\n\n")
    config.write(f"integerisation: {args.integerisation}\n\n")
    config.write(f"bakta_db: {args.bakta_db}\n\n")
    config.write(f"seq_jaccard_threshold: {args.jaccard}\n\n")
    config.write(f"dcj_dist_threshold: {args.dcj}\n\n")
    config.write("prefix: all_plasmids\n\n")
    config.write(f"communities: {args.output_dir}/jaccard/jaccard_communities.txt\n\n")
    config.write(f"dedup_threshold: {args.dedup}\n\n")
    config.write(f"identity_threshold: {args.identity}\n\n")
    config.write(f"timelimit: {timelimit}\n")
    
#integerisation
if args.integerisation == "anno":
    subprocess.run(f"snakemake --snakefile {plingpath}/jac_network_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {forceall}", shell=True, check=True, capture_output=True)
    subprocess.run(f"snakemake --snakefile {plingpath}/anno_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {forceall}", shell=True, check=True, capture_output=True)
elif args.integerisation == "align":
    subprocess.run(f"snakemake --snakefile {plingpath}/align_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {forceall}", shell=True, check=True, capture_output=True)    
else:
    raise Exception(f"\"{args.integerisation}\" is not a valid integerisation method!")
    
#ding and clustering    
subprocess.run(f"snakemake --snakefile {plingpath}/dcj_snakemake/Snakefile --configfile {configfile} --cores {args.cores} --use-conda {forceall}", shell=True, check=True, capture_output=True)


#delete intermediary files
if args.storetmp == False:
    shutil.rmtree(f"{args.output_dir}/tmp_files")
    

