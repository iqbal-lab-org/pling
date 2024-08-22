Advanced Usage
==============

Overview
--------

.. code-block:: console
    
	usage: pling [-h] [--version] [--unimog UNIMOG] [--containment_distance CONTAINMENT_DISTANCE] [--dcj DCJ] [--regions] [--topology TOPOLOGY] [--batch_size BATCH_SIZE] [--sourmash] [--sourmash_threshold SOURMASH_THRESHOLD] [--identity IDENTITY]
		     [--min_indel_size MIN_INDEL_SIZE] [--bh_connectivity BH_CONNECTIVITY] [--bh_neighbours_edge_density BH_NEIGHBOURS_EDGE_DENSITY] [--small_subcommunity_size_threshold SMALL_SUBCOMMUNITY_SIZE_THRESHOLD] [--output_type {html,json,both}]
		     [--plasmid_metadata PLASMID_METADATA] [--ilp_solver {GLPK,gurobi}] [--timelimit TIMELIMIT] [--resources RESOURCES] [--cores CORES] [--profile PROFILE] [--forceall]
		     genomes_list output_dir {align,skip}

	positional arguments:
	  genomes_list          Path to list of fasta file paths.
	  output_dir            Path to output directory.
	  {align,skip}          Integerisation method: "align" for alignment, "skip" to skip integerisation altogether. Make sure to input a unimog file if skipping integerisation.

	options:
	  -h, --help            show this help message and exit
	  --version             show program's version number and exit
	  --unimog UNIMOG       Path to unimog file. Required input if skipping integerisation. (default: None)
	  --containment_distance CONTAINMENT_DISTANCE
		                Threshold for initial containment network. (default: 0.5)
	  --dcj DCJ             Threshold for final DCJ-Indel network. (default: 4)
	  --regions             Cluster regions rather than complete genomes. Assumes regions are taken from circular plasmids. (default: False)
	  --topology TOPOLOGY   File stating whether plasmids are circular or linear. Must be a tsv with two columns, one with plasmid IDs under "plasmid" and one with "linear" or "circular" as entries under "topology". Without this file, pling will asume all plasmids are circular. (default: None)
	  --batch_size BATCH_SIZE
		                How many pairs of genomes to run together in one go (for integerisation from alignment and DCJ calculation steps). (default: 200)
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
	  --output_type {html,json,both}
		                Whether to output networks as html visualisations, cytoscape formatted json, or both. (default: html)
	  --plasmid_metadata PLASMID_METADATA
		                Metadata to add beside plasmid ID on the visualisation graph. Must be a tsv with a single column, with data in the same order as in genomes_list. (default: None)
	  --ilp_solver {GLPK,gurobi}
		                ILP solver to use. Default is GLPK, which is slower but is bundled with pling and is free. If using gurobi, make sure you have a valid license and gurobi_cl is in your PATH. (default: GLPK)
	  --timelimit TIMELIMIT
		                Time limit in seconds for ILP solver. (default: None)
	  --resources RESOURCES
		                tsv stating number of threads and memory to use for each rule. (default: None)
	  --cores CORES         Total number of cores/threads. Put the maximum number of threads you request in the resources tsv here. (This argument is passed on to snakemake's --cores argument.) (default: 1)
	  --profile PROFILE     To run on a cluster with corresponding snakemake profile. (default: None)
	  --forceall            Force snakemake to rerun everything. (default: False)


**Network thresholds:** ``--containment-distance`` and ``--dcj`` control the thresholds at which edges are added to the continament and DCJ-Indel networks. If looking for recent transmissions, we recommend trying threshold 0.3 for the containment distance. If you want to try different DCJ-Indel thresholds with the same containment threshold, you can run pling at different DCJ-Indel thresholds with the same output directory -- this will mean that only the DCJ-Indel network will be calculated anew, which will significantly reduce runtime. Changing the containment threshold will prompt pling to rerun the whole workflow though, and so you should change the output directory if you don't want to lose previous results.

**Genome topologies:** By default, pling assumes all genomes are circular, but with ``--topology`` you can include both circular and linear genomes. Give ``--topology`` a file path to a tsv with the following format:

.. code-block:: console

	plasmid	topology
	genome_1	circular
	genome_2	circular
	genome_3	linear
	genome_4	circular
	
You can also run pling on regions of genomes with the ``--regions`` flag.

**Batching:** ``--batch_size``, ``--sourmash`` and ``--sourmash_threshold`` are all associated with the batching step in pling, in which pairs of plasmids are assigned to a batch. Integerisation and DCJ-Indel calculation is then run per batch, as this improves runtime. If your dataset is in the 1000s, batch sizes of 750 tend to work well. If your dataset is in the order of the 100s or even less, the default batch size should work well enough. If you have a very large (>10k) and diverse dataset, you may want to prefilter which pairs of plasmids you calculate containment distances for and integerise, by first estimating containment distances with sourmash and discarding any too divergent pairs of plasmids early in the workflow.

**Integerisation from alignment parameters:** ``--identity`` and ``--min_indel_size`` control how blocks of sequence are selected during integerisation. For a block of sequence to qualify as shared, its sequence similarity must at least be the value of ``--identity``. Blocks of sequence that are not shared between plasmids are assigned an integer if they are greater than the value of ``--min_indel_size``. Any blocks of sequence less than the value of ``--min_indel_size`` are discarded.

**Clustering parameters:** ``--bh_connectivity`` and ``--bh_neighbour_density`` both determine how a hub plasmid is defined. ``--bh_connectivity`` determines to how many plasmids a hub should at least be connected to, while ``--bh_neighbour_density`` determines how interconnected a hub's neighbours should be.

**Visualisation styles:** With ``--output_type`` you can opt to output cytoscpae compatible json files together with or instead of the default html visualisations. The json files can be loaded directly into cytoscape, and used in python together with networkx. There is also a cytoscape style file available, which provides a style similar to the default html visualisations. It can be downloaded from here: https://raw.githubusercontent.com/leoisl/plasnet/main/plasnet/ext/templates/plasnet_style.xml. You can also add additional metadata to node names with ``--plasmid_metadata``, e.g. Inc types.

**ILP solver:** Calculating the DCJ-Indel distances involves solving and integer linear problem (ILP), and pling allows a choice between two ILP solvers for this: GLPK or gurobi. GLPK is free and bundled with pling, but a bit slower. Gurobi is a commercial software, with a free academic license, and you must have a valid license and installed gurobi beforehand to run pling with it. Both solvers output the same final result. Generally calculating DCJ-Indel is an NP-hard problem, which means in its worst case calculation will take a very long time. The ``--timelimit`` variable sets a time limit for how long the ILP solver takes with a pair, before giving up and outputting the most optimal result it has at that point. However our experience is that when running with integers from alignment, the DCJ-Indel calculation is very quick.

**Snakemake arguments:** Arguments ``--cores``, ``--profile`` and ``--forceall`` are passed as are directly to snakemake. Please refer to snakemake documentation (https://snakemake.readthedocs.io/en/v7.0.0/) for further information. Through ``--resources`` you can pass a path to a ``resources.tsv`` file, which will define number of threads and memory allocated for each rule in pling's snakemake workflows. The format should be the same as the file found under ``pling/resources.tsv``. If you use more than one thread in any rule, remember to set ``--cores`` to the maximum number of threads you'd like to use.

Regions (e.g. integrons)
------------------------

To run pling on a region of a genome, rather than full genomes, use the ``--regions`` flag. Make sure that your fasta files only contain the sequence of your region of interest. This option is intended for comparing specific regions in genomes that undergo rearrangement, e.g. comparing integrons across different genomes. The DCJ-Indel model is not inherently designed for regions of genomes, so pling uses the following assumptions to integerise regions for DCJ-Indel calculation:

1. the sequence outside of the region of interest is identical across all genomes
2. the region of interest is embedded in a circular genome
3. the orientation of the sequences is not known, and should be chosen such that the longest common block of sequence is oriented the same way in both sequences

In practice, assumptions 1 and 2 result in a false integer being added to the integerisation of a pair of sequences. For example, say you have two sequences *A* and *B*, and *1*, *2*, and *3* represent blocks of sequence from *A* and *B*. If *A* and *B* look as follows:

.. code-block:: console

	A: 1 3 2
	B: 1 2

then the integer representation pling uses to calculate DCJ-Indel is:

.. code-block:: console

	A: 1 3 2 4
	B: 1 2 4

Assumption 3 handles the case where e.g. two sequences are the same, but one is the reverse complement of the other. In this case in one of the sequences the false integer has a negative sign, i.e. is reverse complemented, and it would look like this (where *2* is the false integer):

.. code-block:: console

	A: 1 2
	B: -1 -2

The DCJ-Indel distance in this case is then 0, as the genomes containing *A* and *B* are assumed to be circular. Generally, pling decides orientation based off of the longest shared block of sequence. Say for example *A* and *B* have the form:

.. code-block:: console

	A: 1 3 2
	B: 1 -2

where *2* is the longest block of sequence common to A and B. Then the integer representation used by pling would be:

.. code-block:: console

	A: 1 3 2 4
	B: 1 -2 -4

These assumptions were chosen to hopefully best reflect what type of sequences pling will be used on; however, other assumptions may also have been reasonable for modelling regions. It is worth noting though that modifying the assumptions made for calculating on regions usually only slightly change the DCJ-Indel distances (plus or minus a few DCJ-Indel operations), if at all.

You cannot use both the ``--topology`` and ``--regions`` options at the same time -- ``--regions`` will override ``--topology`` and assume all input sequences are regions embedded in circular genomes.

Skipping integerisation and unimog file format
----------------------------------------------

There is the option to skip pling's integerisation and provide your own, by running pling with the ``skip`` option and providing a file path to a unimog formatted file with ``--unimog``. The unimog file will need to contain an integer sequence for each genome, where an integer always denotes the same block of sequence.
Pling performs essentially the same workflow with this option: construction of a containment network, induction of a DCJ-Indel subnetwork from DCJ-Indel distances, and clustering. The only change is that alignment is only used to calculate the containment distances and network, *not* integer sequences, and DCJ-Indel distances are calculated from the user given integerisation, i.e. the unimog file.

The entries in a unimog file are formatted similar to fasta files - the first line begins with ">" and the entry name, and the second lines contains the integer sequence. In the sequence, integers are seperated by a space, and negative integers represent reverse complements. The end of a sequence is indicated by either ")" or "|", where ")" means the genome is circular, and "|" means the genome is linear. Here's an example of what it can look like:

.. code-block:: console

	>genome_1
	1 2 3 )
	>genome_2
	3 1 -2 4 )
	>genome_3
	1 2 4 |
	>genome_4
	2 4 )
	
In this example, there are four entries, and all except for genome_3 are circular.
If you want to construct a unimog file for regions of genomes, rather than full genomes, then you may consider adding a false integer to all the sequences. See the "Regions" section for a discussion of the motivation behind this. In the example from above that would look like this:

.. code-block:: console

	>genome_1
	1 2 3 5 )
	>genome_2
	3 1 -2 4 5 )
	>genome_3
	1 2 4 5 |
	>genome_4
	2 4 5 )

Note that as genome_3 is linear, the false integer here means the region of interest is right at the start of the genome. If this is not the case, you may want to integerise with:

.. code-block:: console

	>genome_3
	6 1 2 4 5 |

for regions that are in the middle of a genome, or:

.. code-block:: console

	>genome_3
	5 1 2 4 |

for regions that are at the end.

For more information on the unimog file format, please refer to https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual.


Improving runtime
-----------------

If you are finding pling is being slow, there's a couple of things you can try that may improve run times:

1. **Increase number of cores:** By default, pling runs with only one core/thread. Increasing the number of cores with ``cores`` can help significantly.
2. **Increase batch size:** Pling assigns pairs of plasmids to a batch, and then integerisation and DCJ-Indel calculation are run per batch, as this improves runtime. You can increase the number of pairs per batch size via ``batch_size``. On larger datasets this can improve runtimes, as increasing the batch sizes (and therefore decreasing the number of batches overall) decreases the overhead of managing those batches. If your dataset is in the 1000s, batch sizes of 750 tend to work well.
3. **Prefilter pairs of plasmids with sourmash:** If you have a very large (>10k) and diverse dataset, you may want to prefilter which pairs of plasmids you calculate containment distances for and integerise, by first estimating containment distances with sourmash and discarding any too divergent pairs of plasmids early in the workflow. Note that sourmash is a kmer sketching approach, and will sometimes fail to correctly estimate containment for very small plasmids (~5kb and smaller). To run the prefiltering, use the ``sourmash`` flag, and you can modify the containment distance at which plasmid pairs get filtered out with ``sourmash_threshold``. We have found that sourmash generally overestimates containment distance, which is why the default threshold for sourmash is 0.85. Decreasing the threshold will filter out more pairs, and therefore decrease run time, but also increases the risk of wrongly discarding pairs. As a general rule of thumb, a sourmash threshold that is higher than the overall containment distance threshold by 0.25 or 0.35 shouldn't wrongly discard pairs.
4. **Use gurobi as the ILP solver:** This is not worth trying unless you have already tried the previous options and still want a quicker run time. By default, pling calculates DCJ-Indel with the ILP solver GLPK, which is free and open source. Gurobi is an alternative, faster ILP solver, but is a commercial software. Please see installation of optional dependancies for more instructions on how to set it up. This is only worth doing if you are planning to use pling repreatedly on large datasets. After installation, use ``--ilp_solver gurobi`` to run pling with gurobi.

Integerisation from annotation (pling v1)
-----------------------------------------

As of version 2, pling no longer does integerisation from gene annotation, as we found it generally unreliable for a variety of reasons. If you are still interested in getting DCJ-Indel distances and clustering from genes rather than pairwise synteny blocks, you can still use version 1. You can download that version from https://github.com/iqbal-lab-org/pling/releases/tag/v1.1.1 , and you will need the following dependancies:

- Python >=3.8, <3.12
- conda
- Snakemake >=7.25.4, <=7.32.4
- pandas >=2.0

To run, required input is a list of paths to fasta files ``genomes_list`` and a path to an output directory ``output_dir``. All the genomes must be circular and complete. If ``pling_path`` is the path to the directory to which you downloaded pling, then usage is::

	PYTHONPATH=<pling_path> python <pling_path>/pling/run_pling.py <genomes_list> <output_dir> anno --bakta_db <bakta_db>

As gene annotation is done via Bakta (https://github.com/oschwengers/bakta), the Bakta database must be downloaded beforehand and provided via ``--bakta_db`` to do integerisation from annotation. If a gene is duplicated multiple times across two plasmids for which you are calculating DCJ-Indel, rather than assigning one integer label to all the paralogs, you may want to match together paralogs that are more similar to each other than the other paralogs. This can speed up the DCJ-Indel claculation, and also provide a more realistic distance. We call this process "deduplication" and it can be controlled via the parameters ``--dedup`` and ``--dedup_threshold``. Note that this approach is scarcely tested, and we did not identified appropriate thresholds, so use at your own risk.
