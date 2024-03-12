# Pling!
Pling is a software workflow for plasmid analysis using rearrangement distances, specifically the Double Cut and Join Indel (DCJ-Indel) distance. By intelligently combining containment distance (shared content as fraction of the smaller) and DCJ-indel distance ("how far apart evolutionarily" in a structural sense), and by preventing shared mobile elements from clouding the issue, it infers clusters of related plasmids.

## Dependancies

- Python >=3.8, <=3.12
- Mamba
- Snakemake >=7.25.4, <=7.32.4
- pandas >=2.0
- Optional: gurobi >=10.0.1

## Installation

```
git clone https://github.com/iqbal-lab-org/pling.git
```

## Basic Usage
Required input is a list of paths to fasta files `genomes_list` and a path to an output directory `output_dir`. All the genomes must be circular and complete. If `pling_path` is the path to the directory to which you downloaded pling, then usage is

```
PYTHONPATH=<pling_path> python <pling_path>/pling/run_pling.py <genomes_list> <output_dir> align
```
for integerisation from alignment (recommended), and
```
PYTHONPATH=<pling_path> python <pling_path>/pling/run_pling.py <genomes_list> <output_dir> anno --bakta_db <bakta_db>
```
for integerisation from annotation, in which `bakta_db` is a path to a Bakta database. For details on the difference between integerisation methods, please see below.

## Description and Output

Pling runs via the python script `run_pling.py`, which creates a config file and runs three snakemake workflows in succession. It starts by calculating containment distances and transforming nucleotide sequences into integer sequences (integerisation, details below), then calculates DCJ-Indel distances, and finally uses these to build a network and cluster on it. The outputs are:

- **Containment communities:** Pling defines broad plasmid communities by building a containment network. In this network, nodes represent plasmids, and an edge is drawn if two plasmids fulfil the containment distance threshold. If plasmid *A* is smaller than plasmid *B*, then the containment distance between the two plasmids is the percentage of plasmid *A* that is *not* contained in plasmid *B*. A plasmid community corresponds exactly to a connected component in this network. Distances and plasmid to community assignments can be found in the folder `containment`.
- **Hub plasmids:** Pling identifies "hub plasmids", which are plasmids that are densely connected on the containment network, but their neighbours are not interconnected. In practice, these are usually relatively small plasmids that consist mostly of a large mobile genetic element, which has spread across many diverse, unrelated plasmids. They are listed in `dcj_thresh_4_graph/objects/hub_plasmids.csv`
- **Integer sequences:** Pling outputs the integer sequences used to calculate DCJ-Indel distances. These are outputted in UniMoG-format (see https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual for description), and also a mapping of integer to sequence coordinates, or a mapping of integer to gene name, is provided. The integer sequences are calculated in batches, so there is a unimog and a map file per batch, all found in the folder `unimogs`.
- **DCJ-Indel subcommunities:** Pling identifies plasmid subcommunities by constructing a DCJ-Indel network, and then clustering on this network. The DCJ-Indel network is a subnetwork of the containment network initially built. Edges are kept if a pair of plasmids fulfil the DCJ-Indel distance threshold. Hub plasmids are isolated in the DCJ-Indel network, with no edges connecting to them (even if they fulfil the DCJ-Indel threshold). On this network Pling clusters using asynchronous label propagation community detection algorithm. Plasmid clusters are labelled by containment community and DCJ-Indel subcommunity, e.g. `community_2_subcommunity_13`, and plasmid to subcommunity assignments are found in `dcj_thresh_4_graph/objects/typing.tsv`. The DCJ-Indel distances are in file `all_plasmid_distances.tsv` (note that the DCJ-Indel distance is only calculated between plasmids which fulfil the containment threshold, so not all pairs of plasmids will be in the file).
- **Visualisations:** Alongside the distances and clustering, Pling outputs network visualisations to aid in further analysis. These can be useful to look if you want to spot interesting relationships between plasmids, e.g. plasmid fusions. They include visualisations of the full containment network and each containment network individually, found under `containment/containment_communities/visualisations`. All nodes are coloured black, but the edges are labelled by containment distance. Additionally, under `dcj_thresh_4_graph/communities` are visualisations of each plasmid community, where nodes are coloured by subcommunity assignment and edges are labelled with both containment distance and DCJ-Indel distance. Finally, under `dcj_thresh_4_graph/subcommunities` are visualisations of each plasmid subcommunity, where nodes all have the same colour, and edges are labelled by containment distance and DCJ-Indel distance. These subcommunity visualisations don't include edges that don't fulfil the DCJ-Indel threshold, but the others do.

## Integerisation

An important step in the workflow is transforming the nucleotide sequences into integer sequences (integerisation), as this is necessary to calculate the DCJ-Indel distances. There are two approaches to this:

**Alignment (recommended):** integers correspond to synteny blocks. The integers are assigned to synteny blocks per plasmid pair, using sequence matches found by nucmer. This approach has been extensively tested and will work on diverse data.

**Annotation:** integers correspond to genes (or conserved blocks of genes). Pangenomes are calculate per plasmid containment community (using panaroo), and this is used for integer assignment. As plasmid communities can be very diverse, panaroo can fail to run on them (this can be somewhat mediated by choosing a stricter containment threshold). Furthermore, gene annotation and pangenome calculation carry a runtime cost, and we have encountered runtime issues with the DCJ-Indel distance calculations as well. If you really want to try this approach, the dataset shouldn't too big and you should know beforehand that there is some similarity between all the plasmids in the dataset. This approach is experimental and we cannot guarantee it will work.

## Advanced Usage

```
usage: run_pling.py [-h] [--version] [--containment_distance CONTAINMENT_DISTANCE] [--dcj DCJ] [--batch_size BATCH_SIZE] [--sourmash] [--sourmash_threshold SOURMASH_THRESHOLD] [--identity IDENTITY]
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
  --version             show program's version number and exit
  --containment_distance CONTAINMENT_DISTANCE
                        Threshold for initial containment network. (default: 0.5)
  --dcj DCJ             Threshold for final DCJ-Indel network. (default: 4)
  --batch_size BATCH_SIZE
                        How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps). (default: 50)
  --sourmash            Run sourmash as first filter on which pairs to calculate DCJ on. Recommended for large and very diverse datasets. (default: False)
  --sourmash_threshold SOURMASH_THRESHOLD
                        Threshold for filtering with sourmash. (default: 0.85)
  --identity IDENTITY   Threshold for percentage of shared sequence between blocks (for integerisation from alignment and for containment calculation). (default: 80)
  --min_indel_size MIN_INDEL_SIZE
                        Minimum size for an indel to be treated as a block (for integerisation from alignment). (default: 200)
  --bh_connectivity BH_CONNECTIVITY
                        Minimum number of connections a plasmid need to be considered a hub plasmid. (default: 10)
  --bh_neighbours_edge_density BH_NEIGHBOURS_EDGE_DENSITY
                        Maximum number of edge density between hub plasmid neighbours to label the plasmid as hub. (default: 0.2)
  --small_subcommunity_size_threshold SMALL_SUBCOMMUNITY_SIZE_THRESHOLD
                        Communities with size up to this parameter will be joined to neighbouring larger subcommunities. (default: 4)
  --plasmid_metadata PLASMID_METADATA
                        Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list. (default: None)
  --ilp_solver {GLPK,gurobi}
                        ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH.
                        (default: GLPK)
  --timelimit TIMELIMIT
                        Time limit in seconds for ILP solver. (default: None)
  --resources RESOURCES
                        tsv stating number of threads and memory to use for each rule. (default: None)
  --cores CORES         Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.) (default:
                        1)
  --profile PROFILE     To run on a cluster with corresponding snakemake profile. (default: None)
  --forceall            Force snakemake to rerun everything. (default: False)
  --dedup               Whether or not to deduplicate (for integerisation from annotation). (default: False)
  --dedup_threshold DEDUP_THRESHOLD
                        Threshold for separating paralogs in deduplication step (for integerisation from annotation). (default: 98.5)
  --bakta_db BAKTA_DB   Path to bakta database (required for integerisation from annotation). (default: None)
```

**Network thresholds:** `--containment-distance` and `--dcj` control the thresholds at which edges are added to the continament and DCJ-Indel networks. If looking for recent transmissions, we recommend trying threshold 0.3 for the containment distance. If you want to try different DCJ-Indel thresholds with the same containment threshold, you can run Pling at different DCJ-Indel thresholds with the same output directory -- this will mean that only the DCJ-Indel network will be calculated anew, which will significantly reduce runtime. Changing the containment threshold will prompt Pling to rerun the whole workflow though, and so you should change the output directory if you don't want to lose previous results.

**Batching:** `--batch_size`, `--sourmash` and `--sourmash_threshold` are all associated with the batching step in Pling, in which pairs of plasmids are assigned to a batch. Integerisation and DCJ-Indel calculation is then run per batch, as this improves runtime. If your dataset is in the 1000s, batch sizes of 200 or 250 tend to work well. If your dataset is in the order of the 100s or even less, the default batch size should work well enough. If you have a very large (>10k) and diverse dataset, you may want to prefilter which pairs of plasmids you calculate containment distances for and integerise, by first estimating containment distances with sourmash and discarding any too divergent pairs of plasmids early in the workflow. Note that currently the sourmash step is memory use heavy, so you may need to adjust resources. See 'Snakemake arguments' below for more information.

**Integerisation from alignment parameters:** `--identity` and `--min_indel_size` control how blocks of sequence are selected during integerisation. For a block of sequence to qualify as shared, its sequence similarity must at least be the value of `--identity`. Blocks of sequence that are not shared between plasmids are assigned an integer if they are greater than the value of `--min_indel_size`. Any blocks of sequence less than the value of `--min_indel_size` are discarded.

**Clustering parameters:** `--bh_connectivity` and `--bh_neighbour_density` both determine how a hub plasmid is defined. `--bh_connectivity` determines to how many plasmids a hub should at least be connected to, while `--bh_neighbour_density` determines how interconnected a hub's neighbours should be.

**ILP solver:** Calculating the DCJ-Indel distances involves solving and integer linear problem (ILP), and Pling allows a choice between two ILP solvers for this: GLPK or gurobi. GLPK is free and bundled with Pling, but a bit slower. Gurobi is a commercial software, with a free academic license, and you must have a valid license and gurobi_cl in your PATH beforehand to run Pling with it. Both solvers output the same final result. Generally calculating DCJ-Indel is an NP-hard problem, which means in its worst case calculation will take a very long time. The `--timelimit` variable sets a time limit for how long the ILP solver takes with a pair, before giving up and outputting the most optimal result it has at that point. However our experience is that when running with integers from alignment, the DCJ-Indel calculation is very quick.

**Snakemake arguments:** Arguments `--cores`, `--profile` and `--forceall` are passed as are directly to snakemake. Please refer to snakemake documentation (https://snakemake.readthedocs.io/en/v7.0.0/) for further information. Through `--resources` you can pass a path to a `resources.tsv` file, which will define number of threads and memory allocated for each rule in Pling's snakemake workflows. The format should be the same as the file found under `pling/resources.tsv`. If you use more than one thread in any rule, remember to set `--cores` to the maximum number of threads you'd like to use.

**Integerisation from annotation parameters:** As gene annotation is done via Bakta (https://github.com/oschwengers/bakta), the Bakta database must be downloaded beforehand and provided via `--bakta_db` to do integerisation from annotation. If a gene is duplicated multiple times across two plasmids for which you are calculating DCJ-Indel, rather than assigning one integer label to all the paralogs, you may want to match together paralogs that are more similar to each other than the other paralogs. This can speed up the DCJ-Indel claculation, and also provide a more realistic distance. We call this process "deduplication" and it can be controlled via the parameters `--dedup` and `--dedup_threshold`. Note that this approach is scarcely tested, and we have not yet identified appropriate thresholds, so use at your own risk.

## Citation

Pling is not yet published.

## Used Tools

- DingII: Gitlab - https://gitlab.ub.uni-bielefeld.de/gi/dingiiofficial; DOI - https://doi.org/10.1007/978-3-030-45257-5_1
- Snakemake 7: Homepage - https://snakemake.readthedocs.io/en/v7.0.0/; DOI - https://doi.org/10.12688/f1000research.29032.1
- Mamba: Homepage - https://mamba.readthedocs.io/en/latest/index.html
- GLPK: Homepage - https://www.gnu.org/software/glpk/
- Gurobi: Homepage - https://www.gurobi.com/
- MUMmer 3.0: Homepage - https://mummer.sourceforge.net/; DOI - https://doi.org/10.1186/gb-2004-5-2-r12
- Sourmash: Github: https://github.com/sourmash-bio/sourmash; DOI - https://doi.org/10.21105/joss.00027
- IntervalTree: Github - https://github.com/chaimleib/intervaltree
- Plasnet: Github - https://github.com/leoisl/plasnet
- NetworkX: Homepage - https://networkx.org/
- Numpy: Homepage - https://numpy.org/; DOI - https://doi.org/10.1038/s41586-020-2649-2
- pandas: Homepage - https://pandas.pydata.org/; DOI - https://zenodo.org/doi/10.5281/zenodo.3509134, https://doi.org/10.25080/Majora-92bf1922-00a
- Biopython: Homepage - https://biopython.org/; DOI - https://doi.org/10.1093/bioinformatics/btp163
- Pafpy: Github - https://github.com/mbhall88/pafpy
- Pymummer: Github - https://github.com/sanger-pathogens/pymummer
- Pyfastaq: Github - https://github.com/sanger-pathogens/Fastaq
- Bakta: Github - https://github.com/oschwengers/bakta; DOI - https://doi.org/10.1099/mgen.0.000685
- Panaroo: Github - https://github.com/gtonkinhill/panaroo; DOI - https://doi.org/10.1186/s13059-020-02090-4
- Minimap2: Github - https://github.com/lh3/minimap2; DOI - https://doi.org/10.1093/bioinformatics/bty191
