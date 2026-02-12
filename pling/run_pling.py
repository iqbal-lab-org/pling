"""
Created on Tue Aug 15 14:34:47 2023

@author: daria
"""

import subprocess
import shutil
import os
import pandas as pd
from pathlib import Path
import yaml
import importlib
from pling import __version__
import click
import logging
from plasnet.utils import PathlibPath
import warnings


def get_pling_path():
    plingpath = os.path.realpath(os.path.dirname(__file__))
    return plingpath

def check_gurobi(ilp_solver):
    if ilp_solver == "gurobi":
        spec = importlib.util.find_spec("gurobi")
        if spec is None:
            raise Exception("Missing optional dependency gurobi!")

def make_config_file(args, integerisation):
    #make configfile
    forceall = ""
    if args["forceall"] == True:
        forceall = "--forceall"

    if args["timelimit"]==None:
        timelimit = "None"
    else:
        if args["ilp_solver"] == "gurobi":
            timelimit= args["timelimit"]
        elif args["ilp_solver"] == "GLPK":
            timelimit = "None"
            warnings.warn("GLPK does not support a time limit; time limit parameter has been ignored.")

    if args["plasmid_metadata"]==None:
        metadata = "None"
    else:
        metadata= str(args["plasmid_metadata"])

    output_dir = Path(args["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = output_dir/"tmp_files"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    configfile = f"{output_dir}/tmp_files/config.yaml"

    config_dict = {"genomes_list": str(args["genomes_list"]), "output_dir": str(output_dir), "integerisation": integerisation, "seq_containment_distance": float(args["containment_distance"]), "dcj_dist_threshold": int(args["dcj"]), "prefix": "all_plasmids","communities": f"{output_dir}/containment/containment_communities", "identity_threshold": float(args["identity"]), "bh_connectivity": int(args["bh_connectivity"]), "bh_neighbours_edge_density": float(args["bh_neighbours_edge_density"]), "small_subcommunity_size_threshold": int(args["small_subcommunity_size_threshold"]), "output_type": str(args["output_type"]), "metadata": metadata, "ilp_solver": str(args["ilp_solver"]), "timelimit": timelimit, "batch_size": int(args["batch_size"])}

    if integerisation=="align":
        config_dict["length_threshold"] = int(args["min_indel_size"])
        if args["topology"]==None:
            config_dict["topology"] = "None"
        else:
            config_dict["topology"]= str(args["topology"])
        if args["regions"]:
            config_dict["regions"] = str(args["regions"])

    if integerisation=="skip":
        config_dict["unimog"] = str(args["unimog"])

    if args["sourmash"]:
        config_dict["sourmash"] = str(args["sourmash"])
        config_dict["sourmash_threshold"] = str(args["sourmash_threshold"])

    if "previous_pling" in args.keys():
        config_dict["reclustering_method"] = args["reclustering_method"]
        if len(args["previous_pling"])==1:
            config_dict["previous_pling"] = str(args["previous_pling"])
        else:
            config_dict["previous_pling"] = ",".join([str(path) for path in args["previous_pling"]])
            if config_dict["reclustering_method"]=="nearest_neighbour":
                raise Exception("Nearest neighbour typing does not support merging graphs.")
    
    profile = ""
    if args["profile"]!=None:
        profile = "--profile " + {args["profile"]}
    if args["resources"]!=None:
        resources = pd.read_csv(args["resources"], sep="\t")
    else:
        resources = pd.read_csv(f"{get_pling_path()}/resources.tsv", sep="\t")
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

def batching(snakemake_args):
    try:
        print("Batching...\n")
        subprocess.run(f"snakemake --snakefile {get_pling_path()}/batching/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
        print("Completed batching.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    
def alignment(snakemake_args):
    try:
        print("Aligning, integerising, and building containment network...\n")
        subprocess.run(f"snakemake --snakefile {get_pling_path()}/align_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
        print("Completed integerisation and containment network.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    
def ding_and_cluster(snakemake_args):
    try:
        print("Calculating DCJ-Indel distance and clustering...\n")
        subprocess.run(f"snakemake --snakefile {get_pling_path()}/dcj_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
        print("Completed distance calculations and clustering.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e



@click.group()
@click.version_option(version=__version__)
@click.help_option()
def cli() -> None:
    pass


@click.group(
        help="Clusters genomes using containment and DCJ-Indel distances."
)
@click.help_option()
def cluster() -> None:
    pass

@click.command(
    help="""
\b
Intergerises genomes using alignment, calculates containment and DCJ-Indel distances, and clusters.
\b
First input is a path to list of fasta file paths.
\b
Second input is a path to an output directory.

""",  # noqa: E501
)
@click.argument("genomes_list", type=PathlibPath(exists=True))
@click.argument("output_dir", type=PathlibPath(exists=False))
@click.option("--containment_distance", default=0.5, help="Threshold for initial containment network.")
@click.option("--dcj", default=4, help="Threshold for final DCJ-Indel network.")
@click.option("--regions", is_flag=True, help="Cluster regions rather than complete genomes. Assumes regions are taken from circular plasmids.")
@click.option("--topology", help="File stating whether plasmids are circular or linear. Must be a tsv with two columns, one with plasmid IDs under \"plasmid\" and one with \"linear\" or \"circular\" as entries under \"topology\". Without this file, pling will asume all plasmids are circular.")
@click.option("--batch_size", default = 200, help="How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps).")
@click.option("--sourmash", is_flag=True, help="Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.")
@click.option("--sourmash_threshold", default=0.85, help="Threshold for filtering with sourmash.")
@click.option("--identity", default=80, help="Threshold for percentage of shared sequence between blocks (for integerisation from alignment and for containment calculation).")
@click.option("--min_indel_size", default=200, help="Minimum size for an indel to be treated as a block (for integerisation from alignment).")
@click.option("--bh_connectivity", default=10, help="Minimum number of connections a plasmid need to be considered a hub plasmid.")
@click.option("--bh_neighbours_edge_density", default=0.2, help="Maximum number of edge density between hub plasmid neighbours to label the plasmid as hub.")
@click.option("--small_subcommunity_size_threshold", default=4, help="Communities with size up to this parameter will be joined to neighbouring larger subcommunities.")
@click.option("--output_type", type=click.Choice(["html", "json", "both"]), default="html", help="Whether to output networks as html visualisations, cytoscape formatted json, or both.")
@click.option("--plasmid_metadata", help="Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list.")
@click.option("--ilp_solver", type=click.Choice(["GLPK", "gurobi"]), default="GLPK",
                    help="ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. "
                            "If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH.")
@click.option("--timelimit", help="Time limit in seconds for ILP solver.")
@click.option("--resources", help="tsv stating number of threads and memory to use for each rule.")
@click.option("--cores", default=1, help="Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.)")
@click.option("--profile", help="To run on a cluster with corresponding snakemake profile.")
#@click.option("--storetmp", is_flag=True, help="Don't delete intermediate temporary files.")
@click.option("--forceall", is_flag=True, help="Force snakemake to rerun everything.")
@click.help_option()
def align(
    **args
):
    check_gurobi(args["ilp_solver"])

    configfile, tmp_dir, forceall, profile = make_config_file(args, "align")

    snakemake_args = f"--configfile {configfile} --cores " + str(args["cores"]) + f" --rerun-incomplete --nolock {profile} {forceall}"

    batching(snakemake_args)

    alignment(snakemake_args)
    
    ding_and_cluster(snakemake_args)

    #delete intermediary files
    shutil.rmtree(tmp_dir)



@click.command(
    help="""
\b    
Given a precomputed integerisation, calculates containment and DCJ-Indel distances, and clusters.
\b
First input is a path to list of fasta file paths.
\b
Second input is a path to an output directory.
\b
Third input is a path to a unimog file. Required input when skipping integerisation.

""",  # noqa: E501
)
@click.argument("genomes_list", type=PathlibPath(exists=True))
@click.argument("output_dir", type=PathlibPath(exists=False))
@click.argument("unimog", type=PathlibPath(exists=True))
@click.option("--containment_distance", default=0.5, help="Threshold for initial containment network.")
@click.option("--dcj", default=4, help="Threshold for final DCJ-Indel network.")
@click.option("--batch_size", default = 200, help="How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps).")
@click.option("--sourmash", is_flag=True, help="Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.")
@click.option("--sourmash_threshold", default=0.85, help="Threshold for filtering with sourmash.")
@click.option("--sourmash_only", is_flag=True, help="Run sourmash instead of aligning to get containment distances. Uses the threshold from \"containment_distance\" rather than \"sourmash_threshold\".")
@click.option("--identity", default=80, help="Threshold for percentage of shared sequence between blocks (for integerisation from alignment and for containment calculation).")
@click.option("--bh_connectivity", default=10, help="Minimum number of connections a plasmid need to be considered a hub plasmid.")
@click.option("--bh_neighbours_edge_density", default=0.2, help="Maximum number of edge density between hub plasmid neighbours to label the plasmid as hub.")
@click.option("--small_subcommunity_size_threshold", default=4, help="Communities with size up to this parameter will be joined to neighbouring larger subcommunities.")
@click.option("--output_type", type=click.Choice(["html", "json", "both"]), default="html", help="Whether to output networks as html visualisations, cytoscape formatted json, or both.")
@click.option("--plasmid_metadata", help="Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list.")
@click.option("--ilp_solver", type=click.Choice(["GLPK", "gurobi"]), default="GLPK",
                    help="ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. "
                            "If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH.")
@click.option("--timelimit", help="Time limit in seconds for ILP solver.")
@click.option("--resources", help="tsv stating number of threads and memory to use for each rule.")
@click.option("--cores", default=1, help="Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.)")
@click.option("--profile", help="To run on a cluster with corresponding snakemake profile.")
#@click.option("--storetmp", is_flag=True, help="Don't delete intermediate temporary files.")
@click.option("--forceall", is_flag=True, help="Force snakemake to rerun everything.")
@click.help_option()
def skip(
    **args
):
    check_gurobi(args["ilp_solver"])

    configfile, tmp_dir, forceall, profile = make_config_file(args, "skip")

    snakemake_args = f"--configfile {configfile} --cores " + str(args["cores"]) + f" --rerun-incomplete --nolock {profile} {forceall}"

    batching(snakemake_args)

    try:
        print("Building containment network...\n")
        subprocess.run(f"snakemake --snakefile {get_pling_path()}/jac_network_snakemake/Snakefile {snakemake_args}", shell=True, check=True, capture_output=True)
        print("Completed containment network.\n")
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e

    ding_and_cluster(snakemake_args)

    #delete intermediary files
    shutil.rmtree(tmp_dir)

cluster.add_command(align)
cluster.add_command(skip)
cli.add_command(cluster)

@click.command(
        help="""
\b
Add new genomes to an existing pling clustering, and/or merge multiple previous clusterings into a single clustering.
\b
First input is a path to list of fasta file paths. This list must contain ALL genomes to be clustered, i.e. both genomes to add to previous clustering and genomes from previous clustering.
\b
Second input is a path to an output directory.
\b
Third input is a path to previous pling output directory (multiple are permitted).

""",  # noqa: E501
)
@click.argument("genomes_list", type=PathlibPath(exists=True))
@click.argument("output_dir", type=PathlibPath(exists=False))
@click.argument("previous_pling", required=True, nargs=-1)
@click.option("--reclustering_method", type=click.Choice(["unbiased", "biased", "nearest_neighbour"]), default="unbiased")
@click.option("--containment_distance", default=0.5, help="Threshold for initial containment network.")
@click.option("--dcj", default=4, help="Threshold for final DCJ-Indel network.")
@click.option("--regions", is_flag=True, help="Cluster regions rather than complete genomes. Assumes regions are taken from circular plasmids.")
@click.option("--topology", help="File stating whether plasmids are circular or linear. Must be a tsv with two columns, one with plasmid IDs under \"plasmid\" and one with \"linear\" or \"circular\" as entries under \"topology\". Without this file, pling will asume all plasmids are circular.")
@click.option("--batch_size", default = 200, help="How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps).")
@click.option("--sourmash", is_flag=True, help="Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.")
@click.option("--sourmash_threshold", default=0.85, help="Threshold for filtering with sourmash.")
@click.option("--identity", default=80, help="Threshold for percentage of shared sequence between blocks (for integerisation from alignment and for containment calculation).")
@click.option("--min_indel_size", default=200, help="Minimum size for an indel to be treated as a block (for integerisation from alignment).")
@click.option("--bh_connectivity", default=10, help="Minimum number of connections a plasmid need to be considered a hub plasmid.")
@click.option("--bh_neighbours_edge_density", default=0.2, help="Maximum number of edge density between hub plasmid neighbours to label the plasmid as hub.")
@click.option("--small_subcommunity_size_threshold", default=4, help="Communities with size up to this parameter will be joined to neighbouring larger subcommunities.")
@click.option("--output_type", type=click.Choice(["html", "json", "both"]), default="html", help="Whether to output networks as html visualisations, cytoscape formatted json, or both.")
@click.option("--plasmid_metadata", help="Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list.")
@click.option("--ilp_solver", type=click.Choice(["GLPK", "gurobi"]), default="GLPK",
                    help="ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. "
                            "If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH.")
@click.option("--timelimit", help="Time limit in seconds for ILP solver.")
@click.option("--resources", help="tsv stating number of threads and memory to use for each rule.")
@click.option("--cores", default=1, help="Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.)")
@click.option("--profile", help="To run on a cluster with corresponding snakemake profile.")
#@click.option("--storetmp", is_flag=True, help="Don't delete intermediate temporary files.")
@click.option("--forceall", is_flag=True, help="Force snakemake to rerun everything.")
@click.help_option()
def add(
    **args
):
    check_gurobi(args["ilp_solver"])

    configfile, tmp_dir, forceall, profile = make_config_file(args, "align")

    snakemake_args = f"--configfile {configfile} --cores " + str(args["cores"]) + f" --rerun-incomplete --nolock {profile} {forceall}"

    batching(snakemake_args)

    alignment(snakemake_args)
    
    ding_and_cluster(snakemake_args)

    #delete intermediary files
    shutil.rmtree(tmp_dir)

cli.add_command(add)

@click.command(
        help = "Output DCJ-Indel distances as matrices for each subcommunity."
)
@click.help_option()
def submatrix():
    pass

cli.add_command(submatrix)

def main():
    cli()

if __name__ == "__main__":
    main()