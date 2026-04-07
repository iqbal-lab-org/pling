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
import sys

logFormatter = logging.Formatter("%(asctime)s [%(levelname)s]  %(message)s")
logger = logging.getLogger(__name__)

def get_pling_path():
    plingpath = os.path.realpath(os.path.dirname(__file__))
    return plingpath

def read_log_file(path, config_dict={}):
    with open(path) as f:
        lines = f.readlines()
        args_bool = True
        i = 4
        while args_bool:
            if lines[i]=="\n":
                args_bool = False
            else:
                arg, value = lines[i].replace(" ", "").strip().split(":")
                config_dict[arg] = value
                i = i + 1
    return config_dict

def check_gurobi(ilp_solver):
    if ilp_solver == "gurobi":
        spec = importlib.util.find_spec("gurobi")
        if spec is None:
            logging.error("Missing optional dependency gurobi!")
            raise Exception("Missing optional dependency gurobi!")
        
def check_vis_trees(vis_trees):
    if vis_trees:
        spec = importlib.util.find_spec("ete3")
        if spec is None:
            logging.warning("Missing optional dependency ete3, will skip tree visualisation.")
            return False
    return True
        
def make_snakemake_config(args, config_dict):
    forceall = ""
    if args["forceall"] == True:
        forceall = "--forceall"

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

    return config_dict, profile, forceall

def set_up_logging(log_file, config_dict):
    with open(log_file, "w") as log:
        log.write(f"pling version {__version__}\n")
        log.write("\n ARGUMENTS \n\n")
        yaml.dump(config_dict,log)
        log.write("\n\n LOG \n\n")
    fileHandler = logging.FileHandler(log_file)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

def make_config_file(args, integerisation):
    #make configfile
    output_dir = Path(os.path.abspath(args["output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = output_dir/"tmp_files"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    configfile = f"{output_dir}/tmp_files/config.yaml"

    config_dict = {"genomes_list": str(args["genomes_list"]), "output_dir": str(output_dir), "integerisation": integerisation, 
                   "seq_containment_distance": float(args["containment_distance"]), "dcj_dist_threshold": int(args["dcj"]), 
                   "prefix": "all_plasmids","communities": f"{output_dir}/containment/containment_communities", 
                   "identity_threshold": float(args["identity"]), "bh_connectivity": int(args["bh_connectivity"]), 
                   "bh_neighbours_edge_density": float(args["bh_neighbours_edge_density"]), "small_subcommunity_size_threshold": int(args["small_subcommunity_size_threshold"]), 
                   "output_type": str(args["output_type"]), "ilp_solver": str(args["ilp_solver"]), "batch_size": int(args["batch_size"]), "visualisation":str(args["visualisation"])}

    if integerisation=="align":
        config_dict["length_threshold"] = int(args["min_indel_size"])
        if args["topology"]==None:
            config_dict["topology"] = "None"
        else:
            config_dict["topology"]= str(os.path.abspath(args["topology"]))
        if args["regions"]:
            config_dict["regions"] = str(args["regions"])

    if integerisation=="skip":
        config_dict["unimog"] = str(args["precomputed"])
        config_dict["sourmash_only"] = str(args["sourmash_only"])

    if args["sourmash"]:
        config_dict["sourmash"] = str(args["sourmash"])
        config_dict["sourmash_threshold"] = str(args["sourmash_threshold"])

    if "previous_pling" in args.keys():
        if args["reclustering_method"]=="asyn":
            config_dict["reclustering_method"] = "unbiased"
        else:
            config_dict["reclustering_method"] = args["reclustering_method"]
        config_dict["previous_pling"] = ",".join([str(os.path.abspath(path)) for path in args["previous_pling"]])
        if config_dict["reclustering_method"]=="nearest_neighbour" and len(args["previous_pling"])>1:
            raise Exception("Nearest neighbour typing does not support merging graphs.")
        for path in args["previous_pling"]:
            prev_thresholds = read_log_file(f"{path}/pling.log")
            if prev_thresholds["dcj_dist_threshold"] != config_dict["dcj_dist_threshold"] or prev_thresholds["seq_containment_distance"] != config_dict["seq_containment_distance"]:
                raise Exception(f"{path} was not constructed with the same containment or DCJ-Indel thresholds as given.")
            
    set_up_logging(f"{output_dir}/pling.log", config_dict)

    if args["timelimit"]==None:
        config_dict["timelimit"] = "None"
    else:
        if args["ilp_solver"] == "gurobi":
            config_dict["timelimit"]= args["timelimit"]
        elif args["ilp_solver"] == "GLPK":
            config_dict["timelimit"] = "None"
            logger.warning("GLPK does not support a time limit; time limit parameter has been ignored.")

    if args["plasmid_metadata"]==None:
        config_dict["metadata"] = "None"
    else:
        config_dict["metadata"]= str(args["plasmid_metadata"])

    config_dict, profile, forceall = make_snakemake_config(args, config_dict)

    with open(configfile, 'w') as config:
        yaml.dump(config_dict, config)

    snakemake_args = f"--configfile {configfile} --cores " + str(args["cores"]) + f" --rerun-incomplete --nolock {profile} {forceall} --quiet all"

    return tmp_dir, snakemake_args

def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        logger.exception(e.stderr.decode())
        raise e

def batching(snakemake_args):
    logger.info("Batching...")
    run_command(f"snakemake --snakefile {get_pling_path()}/batching/Snakefile {snakemake_args}")
    logger.info("Completed batching.")
    
def alignment(snakemake_args):
    logger.info("Aligning, integerising, and building containment network...")
    run_command(f"snakemake --snakefile {get_pling_path()}/align_snakemake/Snakefile {snakemake_args}")
    logger.info("Completed integerisation and containment network.")
    
def ding_and_cluster(snakemake_args):
    logger.info("Calculating DCJ-Indel distance and clustering...")
    run_command(f"snakemake --snakefile {get_pling_path()}/dcj_snakemake/Snakefile {snakemake_args}")
    logger.info("Completed distance calculations and clustering.")


def thresholds(func):
    options = [
        click.option("--containment_distance", default=0.5, help="Threshold for initial containment network."),
        click.option("--dcj", default=4, help="Threshold for final DCJ-Indel network."),
        click.option("--identity", default=80, help="Threshold for percentage of shared sequence between blocks (for integerisation from alignment and for containment calculation)."),
        click.option("--min_indel_size", default=200, help="Minimum size for an indel to be treated as a block (for integerisation from alignment)."),
        click.option("--bh_connectivity", default=10, help="Minimum number of connections a plasmid need to be considered a hub plasmid."),
        click.option("--bh_neighbours_edge_density", default=0.2, help="Maximum number of edge density between hub plasmid neighbours to label the plasmid as hub."),
        click.option("--small_subcommunity_size_threshold", default=4, help="Communities with size up to this parameter will be joined to neighbouring larger subcommunities."),
    ]
    for option in reversed(options):
        func = option(func)
    return func

def vis_paras(func):
    options = [
        click.option("--output_type", type=click.Choice(["html", "json", "both"]), default="html", help="Whether to output networks as html visualisations, cytoscape formatted json, or both."),
        click.option("--plasmid_metadata", help="Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list."),
        click.option("--visualisation", type=click.Choice(["none", "all", "subcommunity"]), default="subcommunity", help="Which network visualisations to produce."),
    ]
    for option in reversed(options):
        func = option(func)
    return func

def resource_paras(func):
    options = [
        click.option("--batch_size", default = 200, help="How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps)."),
        click.option("--ilp_solver", type=click.Choice(["GLPK", "gurobi"]), default="GLPK",
                            help="ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. "
                                    "If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH."),
        click.option("--timelimit", help="Time limit in seconds for ILP solver."),
        click.option("--resources", help="tsv stating number of threads and memory to use for each rule."),
        click.option("--cores", default=1, help="Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.)"),
        click.option("--profile", help="To run on a cluster with corresponding snakemake profile."),
        #@click.option("--storetmp", is_flag=True, help="Don't delete intermediate temporary files.")
        click.option("--forceall", is_flag=True, help="Force snakemake to rerun everything."),
    ]
    for option in reversed(options):
        func = option(func)
    return func

@click.group(context_settings={"show_default": True})
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
@thresholds
@click.option("--regions", is_flag=True, help="Cluster regions rather than complete genomes. Assumes regions are taken from circular plasmids.")
@click.option("--topology", help="File stating whether plasmids are circular or linear. Must be a tsv with two columns, one with plasmid IDs under \"plasmid\" and one with \"linear\" or \"circular\" as entries under \"topology\". Without this file, pling will asume all plasmids are circular.")
@click.option("--sourmash", is_flag=True, help="Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.")
@click.option("--sourmash_threshold", default=0.85, help="Threshold for filtering with sourmash.")
@click.option("--reuse_previous", type=PathlibPath(exists=True), help="Recompute communities and subcommunities from previous containment and DCJ-Indel distances. Provide a path to a previous pling output folder. The containment threshold of previous output MUST have been higher than the current one.")
@vis_paras
@resource_paras
@click.help_option()
def align(
    **args
):
    check_gurobi(args["ilp_solver"])

    tmp_dir, snakemake_args = make_config_file(args, "align")

    if not args["reuse_previous"]:
        batching(snakemake_args)
    else:
        prev_thresholds = read_log_file(args["reuse_previous"]/ "pling.log")
        if float(prev_thresholds["seq_containment_distance"]) >= args["containment_distance"]:
            shutil.copy(args["reuse_previous"] / "containment/all_pairs_containment_distance.tsv", args["output_dir"] / "containment") #copy containment distances - should cause snakemake not to rerun everything, just communities
            shutil.copy(args["reuse_previous"] / "all_plasmids_distances.tsv", args["output_dir"]) #copy DCJ-Indel distances - should cause snakemake not to rerun DCJ-Indel calculations
        else:
            logging.error("Previous containment distance threshold is not higher than current containment distance threshold!")
            raise Exception("Previous containment distance threshold is not higher than current containment distance threshold!")

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
Third input is a path to a unimog file, or previous pling output folder. Required input when skipping integerisation. Unless run with \"reuse_previous\" flag, will assume is a unimog file.

""",  # noqa: E501
)
@click.argument("genomes_list", type=PathlibPath(exists=True))
@click.argument("output_dir", type=PathlibPath(exists=False))
@click.argument("precomputed", type=PathlibPath(exists=True))
@thresholds
@click.option("--sourmash", is_flag=True, help="Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.")
@click.option("--sourmash_threshold", default=0.85, help="Threshold for filtering with sourmash.")
@click.option("--sourmash_only", is_flag=True, help="Run sourmash instead of aligning to get containment distances. Uses the threshold from \"containment_distance\" rather than \"sourmash_threshold\".")
@click.option("--reuse_previous", type=PathlibPath(exists=True), help="Recompute communities and subcommunities from previous containment and DCJ-Indel distances. Provide a path to a previous pling output folder. The containment threshold of previous output MUST have been higher than the current one.")
@vis_paras
@resource_paras
@click.help_option()
def skip(
    **args
):
    check_gurobi(args["ilp_solver"])

    tmp_dir, snakemake_args = make_config_file(args, "skip")

    if not args["reuse_previous"]:
        batching(snakemake_args)
    else:
        prev_thresholds = read_log_file(args["reuse_previous"]/ "pling.log")
        if float(prev_thresholds["seq_containment_distance"]) >= args["containment_distance"]:
            shutil.copy(args["reuse_previous"] / "containment/all_pairs_containment_distance.tsv", args["output_dir"] / "containment") #copy containment distances - should cause snakemake not to rerun everything, just communities
            shutil.copy(args["reuse_previous"] / "all_plasmids_distances.tsv", args["output_dir"]) #copy DCJ-Indel distances - should cause snakemake not to rerun DCJ-Indel calculations
        else:
            logging.error("Previous containment distance threshold is not higher than current containment distance threshold!")
            raise Exception("Previous containment distance threshold is not higher than current containment distance threshold!")

    logger.info("Building containment network...")
    run_command(f"snakemake --snakefile {get_pling_path()}/jac_network_snakemake/Snakefile {snakemake_args}")
    logger.info("Completed containment network.")

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
@click.option("--reclustering_method", type=click.Choice(["asyn", "nearest_neighbour"]), default="asyn")
@thresholds
@click.option("--regions", is_flag=True, help="Cluster regions rather than complete genomes. Assumes regions are taken from circular plasmids.")
@click.option("--topology", help="File stating whether plasmids are circular or linear. Must be a tsv with two columns, one with plasmid IDs under \"plasmid\" and one with \"linear\" or \"circular\" as entries under \"topology\". Without this file, pling will asume all plasmids are circular.")
@click.option("--sourmash", is_flag=True, help="Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.")
@click.option("--sourmash_threshold", default=0.85, help="Threshold for filtering with sourmash.")
@vis_paras
@resource_paras
@click.help_option()
def add(
    **args
):
    check_gurobi(args["ilp_solver"])

    tmp_dir, snakemake_args = make_config_file(args, "align")

    batching(snakemake_args)

    alignment(snakemake_args)
    
    ding_and_cluster(snakemake_args)

    #delete intermediary files
    shutil.rmtree(tmp_dir)

cli.add_command(add)

@click.command(
        help = "Output DCJ-Indel distances as matrices for each subcommunity, and calculate neighbour-joining trees."
)
@click.argument("output_dir", type=PathlibPath(exists=True))
@click.option("--submatrices_dir", type=PathlibPath(exists=False), help="Directory to store results in. Defaults to output_dir/submatrices.")
@click.option("--ignore_containment", is_flag=True, help="Calculate missing DCJ-Indel distances due to pairs not meeting the containment threshold previously. Interpret these with caution!")
@click.option("--vis_trees", is_flag=True, help="Plot the DCJ-Indel NJ trees.")
@resource_paras
@click.help_option()
def submatrix(**args):

    output_dir = str(os.path.abspath(args["output_dir"]))

    config_dict = {}
    if args["submatrices_dir"]:
        config_dict["submatrices_dir"] = str(os.path.abspath(args["submatrices_dir"]))
    else:
        config_dict["submatrices_dir"] = f"{output_dir}/submatrices"
    submatrices_dir = config_dict["submatrices_dir"]
    Path(submatrices_dir).mkdir(parents=True, exist_ok=True)

    config_dict["vis_trees"] = check_vis_trees(args["vis_trees"])
    config_dict["ignore_containment"] = args["ignore_containment"]

    set_up_logging(f"{submatrices_dir}/submatrices.log", config_dict)

    if args["timelimit"]==None:
        config_dict["timelimit"] = "None"
    else:
        if args["ilp_solver"] == "gurobi":
            config_dict["timelimit"]= args["timelimit"]
        elif args["ilp_solver"] == "GLPK":
            config_dict["timelimit"] = "None"
            logger.warning("GLPK does not support a time limit; time limit parameter has been ignored.")

    config_dict = read_log_file(f"{output_dir}/pling.log", config_dict)

    config_dict, profile, forceall = make_snakemake_config(args, config_dict)

    Path(config_dict["submatrices_dir"]).mkdir(parents=True, exist_ok=True)
    configfile = config_dict["submatrices_dir"]+"/config.yaml"
    with open(configfile, 'w') as config:
        yaml.dump(config_dict, config)

    command = "snakemake --cores " + str(args["cores"]) + f" --snakefile {get_pling_path()}/submatrix_snakemake/Snakefile --configfile {configfile} --rerun-incomplete --nolock {profile} {forceall} --quiet all"
    logger.info("Creating distance matrices and trees...")
    run_command(command)
    logger.info("Completed!")

    #delete intermediary files
    os.remove(f"{submatrices_dir}/config.yaml")
    if os.path.exists(f"{submatrices_dir}/incomplete"):
        shutil.rmtree(f"{submatrices_dir}/incomplete")
    if os.path.exists(f"{submatrices_dir}/missing"):
        shutil.rmtree(f"{submatrices_dir}/missing")

cli.add_command(submatrix)

def main():
    cli()

if __name__ == "__main__":
    main()