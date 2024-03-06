# Pling!
Pling is a software workflow for plasmid analysis using rearrangement distances. It calculates containment distances, DCJ-Indel distances and uses these distances to cluster into groups of related plasmids.

## Dependancies

- Python
- Mamba
- Snakemake
- pandas

## Installation

```
git clone https://github.com/iqbal-lab-org/pling.git
```

## Basic Usage
Required input is a list of paths to fasta files `genomes_list` and a path to an output directory `output_dir`. If `pling_path` is the path to the directory to which you downloaded pling, then usage is

```
PYTHONPATH=<pling_path> python <pling_path>/pling/run_pling.py <genomes_list> <output_dir> align
```
for integerisation from alignment (recommended), and
```
PYTHONPATH=<pling_path> python <pling_path>/pling/run_pling.py <genomes_list> <output_dir> anno
```
for integerisation from annotation. For details on the difference between integerisation methods, please see below.

## Description and Output

Pling runs via the python script `run_pling.py`, which creates a config file and runs three snakemake workflows in succession. It starts by calculating containment distances and transforming nucleotide sequences into integer sequences (integerisation, details below), then calculates DCJ-Indel distances, and finally uses these to build a network and cluster on it. The outputs are:

- **Containment network:**
- **Integer sequences:**
- **DCJ-Indel network:**
- **Visualisations:**

## Integerisation

An important step in the workflow is transforming the nucleotide sequences into integer sequences (integerisation), as this is necessary to calculate the DCJ-Indel distances. There are two approaches to this:

**Alignment (recommended):** integers correspond to synteny blocks. The integers are assigned to synteny blocks per plasmid pair, using sequence matches found by nucmer. This approach has been extensively tested and will work on diverse data.

**Annotation:** integers correspond to genes (or conserved blocks of genes). Pangenomes are calculate per plasmid containment community (using panaroo), and this is used for integer assignment. As plasmid communities can be very diverse, panaroo can fail to run on them (this can be somewhat mediated by choosing a stricter containment threshold). Furthermore, gene annotation and pangenome calculation carry a runtime cost, and we have encountered runtime issues with the DCJ-Indel distance calculations as well. If you really want to try this approach, the dataset shouldn't too big and you should know beforehand that there is some similarity between all the plasmids in the dataset. This approach is experimental and we cannot guarantee it will work.

## Advanced Usage

```
usage: run_pling.py [-h] [--containment_distance CONTAINMENT_DISTANCE] [--dcj DCJ] [--batch_size BATCH_SIZE] [--sourmash] [--sourmash_threshold SOURMASH_THRESHOLD] [--identity IDENTITY]
                    [--min_indel_size MIN_INDEL_SIZE] [--bh_connectivity BH_CONNECTIVITY] [--bh_neighbours_edge_density BH_NEIGHBOURS_EDGE_DENSITY]
                    [--small_subcommunity_size_threshold SMALL_SUBCOMMUNITY_SIZE_THRESHOLD] [--plasmid_metadata PLASMID_METADATA] [--ilp_solver {GLPK,gurobi}] [--timelimit TIMELIMIT]
                    [--resources RESOURCES] [--cores CORES] [--profile PROFILE] [--forceall] [--dedup] [--dedup_threshold DEDUP_THRESHOLD] [--bakta_db BAKTA_DB]
                    genomes_list output_dir {anno,align}

positional arguments:
  genomes_list          Path to list of fasta file paths.
  output_dir            Path to output directory.
  {anno,align}          Integerisation method: "anno" for annotation and "align" for alignment. Note that recommended method is integerisation from alignment.

optional arguments:
  -h, --help            show this help message and exit
  --containment_distance CONTAINMENT_DISTANCE
                        Threshold for initial containment network.
  --dcj DCJ             Threshold for final DCJ-Indel network.
  --batch_size BATCH_SIZE
                        How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps).
  --sourmash            Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets.
  --sourmash_threshold SOURMASH_THRESHOLD
                        Threshold for filtering with sourmash.
  --identity IDENTITY   Threshold for percentage of shared sequence between blocks (for integerisation from alignment and for containment calculation).
  --min_indel_size MIN_INDEL_SIZE
                        Minimum size for an indel to be treated as a block (for integerisation from alignment).
  --bh_connectivity BH_CONNECTIVITY
                        Minimum number of connections a plasmid need to be considered a blackhole plasmid.
  --bh_neighbours_edge_density BH_NEIGHBOURS_EDGE_DENSITY
                        Maximum number of edge density between blackhole plasmid neighbours to label the plasmid as blackhole.
  --small_subcommunity_size_threshold SMALL_SUBCOMMUNITY_SIZE_THRESHOLD
                        Communities with size up to this parameter will be joined to neighbouring larger subcommunities.
  --plasmid_metadata PLASMID_METADATA
                        Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list.
  --ilp_solver {GLPK,gurobi}
                        ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH.
  --timelimit TIMELIMIT
                        Time limit in seconds for ILP solver.
  --resources RESOURCES
                        tsv stating number of threads and memory to use for each rule.
  --cores CORES         Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.)
  --profile PROFILE     To run on a cluster with corresponding snakemake profile.
  --forceall            Force snakemake to rerun everything.
  --dedup               Whether or not to deduplicate (for integerisation from annotation).
  --dedup_threshold DEDUP_THRESHOLD
                        Threshold for separating paralogs in deduplication step (for integerisation from annotation).
  --bakta_db BAKTA_DB   Path to bakta database (required for integerisation from annotation).

```

**Network thresholds:**

**Batching:**

**Integerisation from alignment parameters:**

**Clustering parameters:**

**ILP solver:**

**Snakemake arguments:** 

**Integerisation from annotation parameters:**

## Citation
